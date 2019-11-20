program wwz
    !----------------------------------------------------------------------------
    !----------------------------------------------------------------------------
    !
    !     WEIGHTED WAVELET Z-TRANSFORM
    !
    !     This is a fortran version of Grant Foster's WWZ1.1.BAS BASIC program,
    !     stripped down for speed and readability improvements.  It is modestly
    !     less flexible than the basic code (esp regarding varying input file
    !     formats), but will be easy to modify to suit your needs.
    !
    !     In the event you have more than 100000 data points, you can resize all
    !     arrays to meet your needs -- the only limitation is your free memory.
    !     (Note that all variables are double precision, though.)
    !
    !     A description of the mathematics can be found in G. Foster, "Wavelets
    !     for Period Analysis of Unevenly Sampled Time Series", Astronomical Journal
    !     112, 1709 (Oct 1996).
    !
    !		-Matthew Templeton, August 15, 2002
    !----------------------------------------------------------------------------

    ! print copyright statement from original program

    write(6, 200)
    200   format(/////, "    Weighted Wavelet Z-transform (WWZ)", /, &
            "  (C) Copyright 1996, 2002 by the American Association", /, &
            "    of Variable Star Observers; all rights reserved."&
            , ////)

    ! main program
    call dataread
    call maketau
    call makefreq
    call getcon
    call writehead
    call wwt
    write(6, *) "Program complete!"
    stop
    ! end of program
end


subroutine dataread
    ! reads in a data file
    ! note: this assumes your data are in two columns: time(1) and mag(2)
    ! it can be easily changed by buffering columns with dummy variables -- just
    ! declare "real*8 dum" and read in null columns with that variable.

    ! this subroutine also returns the number of data points (numdat)
    ! the average (dave), the variance (dvar), and the standard deviation (dsig),
    ! along with the name of the input file.

    implicit none
    integer i, numdat
    character*80 filen, fileo
    real*8 dt, dx
    real*8 dave, dvar, dsig
    common/stardat/dt(100000), dx(100000), dave, dvar, dsig, numdat
    common/outname/filen

    print*, 'input data filename'
    read*, filen
    open(unit = 1, file = filen, status = 'old')
    numdat = 0
    do i = 1, 100000
        read(1, *, end = 999)dt(i), dx(i)
        dave = dave + dx(i)
        dvar = dvar + (dx(i) * dx(i))
        numdat = numdat + 1
    enddo
    ! when eof(1) is reached, the loop jumps here
    999   continue
    close(1)
    dave = dave / dfloat(numdat)
    dvar = (dvar / dfloat(numdat)) - (dave * dave)
    dsig = dsqrt(dvar * dfloat(numdat) / dfloat(numdat - 1))

    ! open your output file (there is only one output file in this version)
    print*, 'output wwz filename (e.g. name.wwz)'
    read*, fileo
    open(unit = 2, file = fileo, status = "unknown")
    open(unit = 3, file = 'wwzper.dat', status = "unknown")
    return
end


subroutine maketau
    ! make your array of time lags, tau, here.  In BASIC, you could just step from
    ! whatever taus you choose, like "FOR DTAU = DTAU1 TO DTAU2 STEP DT", but FORTRAN
    ! only lets you loop over integers.  So we determine the taus ahead of time.

    implicit none
    real*8 tau, dt, dx
    real*8 dave, dvar, dsig
    real*8 dtspan, dtstep, dtaulo, dtauhi
    real*8 taucheck
    integer numdat, i, j, ntau

    common/stardat/dt(100000), dx(100000), dave, dvar, dsig, numdat
    common/taudat/tau(100000), ntau

    dtspan = dt(numdat) - dt(1)
    dtstep = dtspan / 50.d0
    call round(dtstep)
    dtaulo = dt(1)
    dtauhi = dt(numdat)
    dtaulo = dtstep * dfloat(idint((dtaulo / dtstep) + 0.5d0))
    dtauhi = dtstep * dfloat(idint((dtauhi / dtstep) + 0.5d0))
    tau(1) = dtaulo
    ntau = 1
    do i = 2, 100000
        taucheck = tau(1) + dfloat(i - 1) * dtstep
        if(taucheck.gt.dtauhi) goto 998
        tau(i) = taucheck
        ntau = ntau + 1
    enddo
    998   continue
    return
end


subroutine round(darg)
    ! rounds the taus... from G. Foster's code.
    implicit none
    real*8 dex, nex, darg
    dex = dlog10(darg)
    nex = dint(dex)
    darg = darg / (10.**nex)
    if (darg.ge.5.d0) then
        darg = 5.d0
    else
        if (darg.ge.2.d0) then
            darg = 2.d0
        else
            darg = 1.d0
        endif
    endif
    darg = darg * (10**nex)
    return
end


subroutine makefreq
    ! query the user for a frequency range, and make an array of frequencies
    implicit none
    real*8 flo, fhi, df, freq
    integer nfreq, i
    character*1 choice
    common/freqdat/freq(100000), nfreq

    700   continue
    print*, 'low frequency (cyc/d)'
    read*, flo
    print*, 'high frequency'
    read*, fhi
    print*, 'delta f'
    read*, df
    nfreq = nint((fhi - flo) / df) + 1

    ! depending upon the user inputs, you may wind up with a lot of frequencies.
    ! this will query in case nfreq>1000 (nfreq=1000 is a reasonable number)
    if(nfreq.gt.1000) then
        701    print*, 'you have more than 1000 frequencies.  OK? (y/n)'
        read*, choice
        if(choice.eq.'N'.or.choice.eq.'n') goto 700
        if(choice.ne.'Y'.and.choice.ne.'y') goto 701
    endif

    freq(1) = flo
    do i = 2, nfreq
        freq(i) = freq(1) + dfloat(i - 1) * df
    enddo

    return
end


subroutine getcon
    ! query for a decay constant, C (as in exp(-c*omega^2)), for the wavelet window
    implicit none
    real*8 dcon
    common/condat/dcon
    dcon = 0.001
    print*, 'input decay constant c'
    read*, dcon
    return
end


subroutine matinv
    ! invert the matrix of the wwz equations...
    implicit none
    real*8 dmat
    real*8 dsol(0:2, 0:2)
    real*8 dfac
    integer ndim, i, j, k, ni, nj
    common/matdat/dmat(0:2, 0:2)

    ndim = 2

    do i = 0, 2
        do j = 0, 2
            dsol(i, j) = 0.d0
        enddo
        dsol(i, i) = 1.d0
    enddo

    do i = 0, ndim
        if(dmat(i, i).eq.0.d0) then
            if (i.eq.ndim) return
            do j = i + 1, ndim
                if(dmat(j, i).ne.0.d0) then
                    do k = 0, ndim
                        dmat(i, k) = dmat(i, k) + dmat(j, k)
                        dsol(i, j) = dsol(i, j) + dsol(j, k)
                    enddo
                endif
            enddo
        endif
        dfac = dmat(i, i)
        do j = 0, ndim
            dmat(i, j) = dmat(i, j) / dfac
            dsol(i, j) = dsol(i, j) / dfac
        enddo
        do j = 0, ndim
            if(j.ne.i) then
                dfac = dmat(j, i)
                do k = 0, ndim
                    dmat(j, k) = dmat(j, k) - (dmat(i, k) * dfac)
                    dsol(j, k) = dsol(j, k) - (dsol(i, k) * dfac)
                enddo
            endif
        enddo
    enddo
    do i = 0, ndim
        do j = 0, ndim
            dmat(i, j) = dsol(i, j)
        enddo
    enddo
    return
end


subroutine writehead
    ! Write a similar header for this output file as is done for the formatted
    ! output of the BASIC code.
    implicit none
    character*80 filen
    real*8 tau, freq
    real*8 dt, dx
    real*8 dave, dvar, dsig
    integer numdat, ntau, nfreq
    common/outname/filen
    common/taudat/tau(100000), ntau
    common/freqdat/freq(100000), nfreq
    common/stardat/dt(100000), dx(100000), dave, dvar, dsig, numdat

    write(2, 201) filen, numdat, dave, dsig, dvar
    write(2, 202) tau(1), tau(ntau), tau(2) - tau(1)
    write(2, 203) freq(1), freq(nfreq), freq(2) - freq(1)
    write(2, 204)

    201   format("File=", a13, " NUM=", i6, " AVE=", f11.4, " SDV=", f11.4, " VAR=", f11.4)
    202   format("    From JD ", f13.4, " to JD ", f13.4, " step", f11.4)
    203   format("   From fre", f11.7, " to fre", f11.7, " step", f11.7)
    204   format(8h     tau, 13h         freq, 16h           WWZ, 12h         Amp, 14h        m(ave), 12h        Neff)
    return
end


subroutine wwt
    implicit none
    real*8 dt, dx, tau, freq, twopi, dcon
    real*8 dvec(0:2), dmat, dcoef(0:2)
    real*8 dave, dvar, dsig
    integer ndim
    integer itau1, itau2, ifreq1, ifreq2, nfreq, ntau
    integer i, j, itau, ifreq, idat, nstart, numdat
    real*8 domega, dweight2, dz, dweight
    real*8 dcc, dcw, dss, dsw, dxw, dvarw
    real*8 dtau
    real*8 dpower, dpowz, damp, dneff, davew
    real*8 dfre
    ! dnefff should probably be dneff, but we'll keep it for strict consistency...
    real*8 dnefff
    integer n1, n2

    real*8 dmz, dmfre, dmper, dmpow, dmamp, dmcon, dmneff

    common/stardat/dt(100000), dx(100000), dave, dvar, dsig, numdat
    common/taudat/tau(100000), ntau
    common/freqdat/freq(100000), nfreq
    common/constants/twopi
    common/condat/dcon
    common/matdat/dmat(0:2, 0:2)

    twopi = 2.d0 * dacos(-1.d0)
    ndim = 2
    itau1 = 1
    itau2 = ntau
    ifreq1 = 1
    ifreq2 = nfreq
    nstart = 1

    do itau = itau1, itau2
        nstart = 1
        dtau = tau(itau)

        dmz = 0.d0
        do ifreq = ifreq1, ifreq2
            dfre = freq(ifreq)
            domega = dfre * twopi
            do i = 0, ndim
                dvec(i) = 0.d0
                do j = 0, ndim
                    dmat(i, j) = 0.d0
                enddo
            enddo
            dweight2 = 0.d0

            do idat = nstart, numdat
                dz = domega * (dt(idat) - dtau)
                dweight = dexp(-1.d0 * dcon * dz * dz)
                if (dweight.gt.1.d-9) then
                    dcc = dcos(dz)
                    dcw = dweight * dcc
                    dss = dsin(dz)
                    dsw = dweight * dss
                    dmat(0, 0) = dmat(0, 0) + dweight
                    dweight2 = dweight2 + (dweight * dweight)
                    dmat(0, 1) = dmat(0, 1) + dcw
                    dmat(0, 2) = dmat(0, 2) + dsw
                    dmat(1, 1) = dmat(1, 1) + (dcw * dcc)
                    dmat(1, 2) = dmat(1, 2) + (dcw * dss)
                    dmat(2, 2) = dmat(2, 2) + (dsw * dss)
                    dxw = dweight * dx(idat)
                    dvec(0) = dvec(0) + dxw
                    dvarw = dvarw + (dxw * dx(idat))
                    dvec(1) = dvec(1) + (dcw * dx(idat))
                    dvec(2) = dvec(2) + (dsw * dx(idat))
                else
                    if (dz.gt.0.d0) then
                        goto 902
                    else
                        nstart = idat + 1
                    endif
                endif
            enddo

            ! go here once dz>0 (In BASIC version, this is the "OUTLOOP" section of SUB WWZ)
            902     continue
            dpower = 0.d0
            damp = 0.d0
            do n1 = 0, ndim
                dcoef(n1) = 0.d0
            enddo
            if(dweight2.gt.0.d0) then
                dneff = (dmat(0, 0) * dmat(0, 0)) / dweight2
            else
                dneff = 0.d0
            endif
            if(dneff.gt.3.d0) then
                do n1 = 0, ndim
                    dvec(n1) = dvec(n1) / dmat(0, 0)
                    do n2 = 1, ndim
                        dmat(n1, n2) = dmat(n1, n2) / dmat(0, 0)
                    enddo
                enddo
                if (dmat(0, 0).gt.0.d0) then
                    dvarw = dvarw / dmat(0, 0)
                else
                    dvarw = 0.d0
                endif
                dmat(0, 0) = 1.d0
                davew = dvec(0)
                dvarw = dvarw - (davew * davew)
                if (dvarw.le.0.d0) dvarw = 1.d-12
                do n1 = 1, ndim
                    do n2 = 0, n1 - 1
                        dmat(n1, n2) = dmat(n2, n1)
                    enddo
                enddo

                call matinv

                do n1 = 0, ndim
                    do n2 = 0, ndim
                        dcoef(n1) = dcoef(n1) + dmat(n1, n2) * dvec(n2)
                    enddo
                    dpower = dpower + (dcoef(n1) * dvec(n1))
                enddo
                dpower = dpower - (davew * davew)
                dpowz = (dneff - 3.d0) * dpower / (dvarw - dpower) / 2.d0
                dpower = (dneff - 1.d0) * dpower / dvarw / 2.d0
                damp = dsqrt((dcoef(1) * dcoef(1)) + (dcoef(2) * dcoef(2)))
            else
                dpowz = 0.d0
                dpower = 0.d0
                damp = 0.d0
                if (dneff.lt.1.d-9) dnefff = 0.d0
            endif
            if (damp.lt.1.d-9) damp = 0.d0
            if (dpower.lt.1.d-9) dpower = 0.d0
            if (dpowz.lt.1.d-9) dpowz = 0.d0

            ! now, write everything out -- one write per frequency per tau...
            write(2, 205) dtau, dfre, dpowz, damp, dcoef(0), dneff
            205     format(f12.4, 2x, f10.7, 4(2x, f11.4))
            if(dpowz.gt.dmz) then
                dmz = dpowz
                dmfre = dfre
                dmper = 1.d0 / dfre
                dmpow = dpower
                dmamp = damp
                dmcon = dcoef(0)
                dmneff = dneff
            endif

        enddo
        write(3, 206) dtau, dmper, dmamp, dmcon, dmfre, dmz, dmneff
        206    format(f13.4, f11.4, f14.4, f11.4, f11.7, 2(f11.4))

    enddo
    ! ...and close the output file when you're done.
    close(2)
    close(3)
    return
end