program elektron2dsi
    use, intrinsic :: iso_c_binding
    use fourier
    use fun
    use fourier

    implicit none

    !type(c_ptr), intent(in), value :: ptr
    real(c_double) eff_tmp, eff_tmp_k, eff_tmp_b, eff_tmp_k_b, eff(2), &
        sum_eff, loss_on_the_way_plus, loss_on_the_way_plus_k, loss_on_the_way_minus, &
     int_abs2_a_plus_at_z0, int_abs2_a_plus_at_z0_k, int_abs2_a_plus_at_zl, int_abs2_a_plus_at_zl_k, int_abs2_a_plus_at_zl_on_mir, &
        int_abs2_a_plus_at_zl_out_mir, int_abs2_a_plus_at_zl_out_mir_k, int_abs2_a_minus_at_z0_out_mir_k, &
            int_abs2_a_minus_at_z0_k, int_abs2_a_minus_at_z0_on_mir, int_abs2_a_minus_at_z0_out_mir, int_abs2_a_minus_at_zl_on_mir, loss_on_the_way_minus_k, &
        int_abs2_a_minus_at_zl_on_mir_k, int_obrezki_z0, int_obrezki_zl, int_abs2_a_minus_at_z0, int_abs2_a_minus_at_zl, obrezki
    integer(c_int) iz, percent, rec_length, nr, first_it, ith
    character*6 str
    integer i

    !chk_arr = check_alloc()
    !if (chk_arr /= .true.) stop 'allocation error!'

    !read out the parameters, calculate a0(x), f(x), k2, z, x, etc.
    call init()

    !what does not depend on j
    dlt = im1*gamma*k2 + sigma
    ex = cdexp(-dlt*h)

    rec_length = c_double_complex*2*nx/4
    !open files
    open (221, file='$$$z0.bin', access='direct', recl=rec_length, err=101)
    open (222, file='$$$zl.bin', access='direct', recl=rec_length, err=101)
    open (3, file='eff.dat', err=101)
    open (35, file='eff_b.dat', err=101)
    open (53, file='eff_new.dat', err=101)
    open (33, file='eff_k.dat', err=101)
    open (335, file='eff_k_b.dat', err=101)
    open (533, file='eff_new_k.dat', err=101)
    open (788, file='it.dat', err=101)
    open (777, file='for_graphics.dat', err=101)

    call cpu_time(start_time)

    first_it = it_made + 1
    do_iter: do it_made = first_it, it_doiter
        percent = 0
        write (*, '(a11,\,i5,\,a2,i5,a,\)') 'iteration #', it_made, ': ', percent, '%'

        !initial conditions for the phase and its distunction at z=0
        call calc_theta(th0, dthdz)

        eff(1) = 1.0d0/nth*sum(dthdz(:, 1) - delta) !(xcp1)
        eff(2) = 1.0d0/nth*sum(dthdz(:, 2) - delta) !(xcp2)

        !initial current (z=0)
        jk0 = fk1*jf(th0(:, 1)) + fk2*jf(th0(:, 2))

        !ak0
        !ak0 = fs(a0)
        call makeak(ak0, a0)

        !setting the file index
        inquire (222, nextrec=nr)
        !write a(z=0,x)
        write (221, rec=nr) a0
        !print *, 'nr = ', nr

        if (((intrvl > 0) .and. (mod(it_made, intrvl) == 0)) &
            .or. (it_made == first_it) &
            .or. (it_made == it_doiter)) then
            a_amp_z0 = cdabs(a0)
            a_spec_amp_z0 = cdabs(fs(a0))
            write (str, 105) it_made
105         format(i0.6)
            str = trim(str)
            open (1, file='a_'//str//'.bin', form='binary', err=101)
            open (2, file='ak_'//str//'.bin', form='binary', err=101)
            open (88, file='jk_'//str//'.bin', form='binary', err=101)
            write (1, err=103) a0z0
            write (2, err=103) ak0
            write (88, err=103) jk0
            if (amp_only == .false.) then
                open (8, file='eff_'//str//'.bin', form='binary', err=101)
                write (8, err=103) eff
            end if
            if (thout == .true.) then
                open (10, file='th_'//str//'.dat', err=101)
                write (10, '(e17.8,\)', err=103) 0.0
                do ith = 1, nth
                    write (10, '(e17.8,\)', err=103) th0(ith, 1)
                end do
                do ith = 1, nth
                    write (10, '(e17.8,\)', err=103) th0(ith, 2)
                end do
                write (10, '(/,\)')
            end if
        end if

        !to calculate the efficiency we calculate the sum of abs(a(z=0))**2 by x
        int_abs2_a_plus_at_z0 = sum(cdabs(a0)*cdabs(a0))*hx
        int_abs2_a_plus_at_z0_k = 0.5d0*sum(cdabs(ak0)*cdabs(ak0))

        !sum of field amplitudes at z=0 (also for effectiveness)
        sum_abs2_a_plus_by_z = cdabs(a0)*cdabs(a0)
        sum_abs2_a_plus_by_z_k = 0.5d0*cdabs(ak0)*cdabs(ak0)

        do_z: do iz = 1, nz - 1
            rhs0 = rhs(ak0, th0)
            th1 = th0 + dthdz*h + h/2.0d0*rhs0*h !theta predictor
            jk1 = fk1*jf(th1(:, 1)) + fk2*jf(th1(:, 2)) !current predictor

            !predictor a (interpolation)
            ak1(1) = ak0(1) + h/2.0d0*(jk0(1) + jk1(1))
            ak1(2:nk) = ak0(2:nk)*ex(2:nk) + &
                        c3*(jk0(2:nk) + jk1(2:nk)*(-1.0d0 + dlt(2:nk)*h) + &
                            ex(2:nk)*(jk1(2:nk) - jk0(2:nk)*(1.0d0 + dlt(2:nk)*h)))/dlt(2:nk)/dlt(2:nk)/h

            !predictor a (keystone)
            !atmp = (ak0 + c3 * h / 2.0d0 * jk0) * cdexp(-dlt * h) !part a
            !ak1 = atmp + c3 * h / 2.0d0 * jk1 !predictor a

            !a1 = ifs(ak1) !back to reality
            call makea(a1, ak1)

            !corrector theta - 1
            th1 = th0 + dthdz*h + h/6.0d0*rhs0*h &
                  + h/3.0d0*rhs(ak1, th1)*h

            !theta(:, :, iz + 1) = th1 !write down the electron trajectory

            jk1 = fk1*jf(th1(:, 1)) + fk2*jf(th1(:, 2)) !current corrector

            !what depends on the current
            !jkd = jk1 - jk0

            !offset a (interpolation)
            ak1(1) = ak0(1) + h/2.0d0*(jk0(1) + jk1(1))
            ak1(2:nk) = ak0(2:nk)*ex(2:nk) + &
                        c3*(jk0(2:nk) + jk1(2:nk)*(-1.0d0 + dlt(2:nk)*h) + &
                            ex(2:nk)*(jk1(2:nk) - jk0(2:nk)*(1.0d0 + dlt(2:nk)*h)))/dlt(2:nk)/dlt(2:nk)/h

            !corrector a (keystone)
            !atmp = (ak0 + c3 * h / 2.0d0 * jk0) * cdexp(-dlt * h) !part a
            !ak1 = atmp + c3 * h / 2.0d0 * jk1 !corrector a

            dthdz = dthdz + h/2.0d0*(rhs0 + rhs(ak1, th1))
            !corrector theta - 1 end

            !calculation of efficiency at point z = iz*h
            !doix=1,nx
            ! eff(ix) = 1.0d0 / nth * sum(dthdz(:,ix) - delta) !counting efficiency
            !end do
            eff(1) = 1.0d0/nth*sum(dthdz(:, 1) - delta) !counting efficiency (xcp1)
            eff(2) = 1.0d0/nth*sum(dthdz(:, 2) - delta) !counting efficiency (xcp2)

            !back to reality
            !a1 = ifs(ak1)
            call makea(a1, ak1)

            !sum of field amplitudes at z=iz*h
            sum_abs2_a_plus_by_z = sum_abs2_a_plus_by_z + cdabs(a1)*cdabs(a1)
            sum_abs2_a_plus_by_z_k = sum_abs2_a_plus_by_z_k + 0.5d0*cdabs(ak1)*cdabs(ak1)

            if (((intrvl > 0) .and. (mod(it_made, intrvl) == 0)) &
                .or. (it_made == first_it) &
                .or. (it_made == it_doiter)) then
                write (1, err=103) a1
                write (2, err=103) ak1
                write (88, err=103) jk1
                if (amp_only == .false.) then
                    write (8, err=103) eff
                end if
                if (thout == .true.) then
                    write (10, '(e17.8,\)', err=103) iz*h
                    do ith = 1, nth
                        write (10, '(e17.8,\)', err=103) th1(ith, 1)
                    end do
                    do ith = 1, nth
                        write (10, '(e17.8,\)', err=103) th1(ith, 2)
                    end do
                    write (10, '(/,\)')
                end if
            end if

            th0 = th1 !for the next z step
            jk0 = jk1 !for the next step in z
            ak0 = ak1 !for the next step in z
            !dthdz0 = dthdz1 !for the next z step

            percent = int(real(iz)/real(nz - 2)*100 + 0.5)
            write (*, '(\,a6,i5,a)') '\b\b\b\b\b\b'c, percent, '%'
        end do do_z

        !calculation of the sum of efficiency for x at the point z=lz
        sum_eff = (eff(1) + eff(2))/2.0d0

        !to calculate the efficiency we calculate the sum of abs(a(z=lz))**2 by x
        int_abs2_a_plus_at_zl = sum(cdabs(a1)*cdabs(a1))*hx
        int_abs2_a_plus_at_zl_k = 0.5d0*sum(cdabs(ak1)*cdabs(ak1))

        int_abs2_a_plus_at_zl_on_mir = sum(cdabs(a1(ig0:ig1 - 1))*cdabs(a1(ig0:ig1 - 1)))*hx
        int_abs2_a_plus_at_zl_out_mir = (sum(cdabs(a1(1:ig0 - 1))*cdabs(a1(1:ig0 - 1))) &
                                         + sum(cdabs(a1(ig1 + 1:nx))*cdabs(a1(ig1 + 1:nx))))*hx

        aktmp = a1*g
        call sint(aktmp)
        !subtract the uncut modes of the cut reflection
        int_abs2_a_plus_at_zl_out_mir_k = int_abs2_a_plus_at_zl_k - 0.5d0*sum(cdabs(aktmp)*cdabs(aktmp))

        aktmp(1:nk) = dcmplx(0)
        int_obrezki_zl = 0.5d0*dmysum(cdabs(aktmp)*cdabs(aktmp)) ! energy in cut to zl

        !double sum of field amplitudes in z and in x
        loss_on_the_way_plus = 2.0d0*sigma*sum(sum_abs2_a_plus_by_z)*hx*h
        loss_on_the_way_plus_k = 2.0d0*sigma*sum(sum_abs2_a_plus_by_z_k)*h

        !write a(z=l,x)
        write (222, rec=nr) a1
        write (788, *) it_made

        if (it_made == it_doiter) then
            a_amp_zl = cdabs(a1)
            a_spec_amp_zl = cdabs(fs(a1))
        end if

        !initial conditions for the next iteration
        !if (it_made < it_doiter) then
        if (recount == .false.) then
            a0 = r1*a1*g*r0
        else
            atmp = r1*a1*g!! a minus on aperture z \u003d lz (real)

            !fourier
            call makeak(akzl, atmp)

            call makea(atmp, akzl)

            int_abs2_a_minus_at_zl_on_mir_k = 0.5d0*dmysum(cdabs(akzl)*cdabs(akzl)) !for efficiency (what is already reflected from the mirror)
            int_abs2_a_minus_at_zl_on_mir = sum(cdabs(atmp)*cdabs(atmp))*hx

            !------------------------------------------------- -------------------------------------------------- -------------------------------------------------- --

            akz0 = akzl*cdexp(-dlt*lz) !what is returned to z0

            int_abs2_a_minus_at_z0_k = 0.5d0*dmysum(cdabs(akz0)*cdabs(akz0)) !a minus on aperture z=0 (real)

            call make_a0z0_ak1_atmp(a0z0, a0z0cut, ak0, akz0)

            int_abs2_a_minus_at_z0 = sum(cdabs(a0z0)*cdabs(a0z0))*hx

            !akz0 fashions of what is going to the beginning
            !ak0 cut mode of reflection cut

            loss_on_the_way_minus_k = int_abs2_a_minus_at_zl_on_mir_k - int_abs2_a_minus_at_z0_k

            aktmp = a0z0*g
            call sint(aktmp) ! uncut reflection cut modes
            int_abs2_a_minus_at_z0_out_mir_k = int_abs2_a_minus_at_z0_k - 0.5d0*sum(cdabs(aktmp)*cdabs(aktmp))

            aktmp(1:nk) = dcmplx(0)
            int_obrezki_z0 = 0.5d0*dmysum(cdabs(aktmp)*cdabs(aktmp)) ! energy in cut to z0

            int_abs2_a_minus_at_z0_on_mir = sum(cdabs(a0z0(ig0:ig1 - 1))*cdabs(a0z0(ig0:ig1 - 1)))*hx
                    int_abs2_a_minus_at_z0_out_mir = (sum(cdabs(a0z0(1:ig0-1))*cdabs(a0z0(1:ig0-1))) + sum(cdabs(a0z0(ig1+1:nx))*cdabs(a0z0(ig1+ 1:nx)))) * hx
            loss_on_the_way_minus = int_abs2_a_minus_at_zl_on_mir - (int_abs2_a_minus_at_z0_on_mir + int_abs2_a_minus_at_z0_out_mir)

            !reflection and cut at the first mirror
            a0 = a0z0cut
            open (388, file='a0.bin', form='binary', err=101)
            write (388) a0
            close (388)
        end if
        !endif

        !recording efficiency at this iteration
        eff_tmp = (loss_on_the_way_plus + int_abs2_a_plus_at_zl - int_abs2_a_plus_at_z0)/lx - 4.0d0*c3*sum_eff
        eff_tmp_k = loss_on_the_way_plus_k + int_abs2_a_plus_at_zl_k - int_abs2_a_plus_at_z0_k - 4.0d0*c3*sum_eff

        obrezki = int_obrezki_z0 + int_obrezki_zl

        write (3, 104, err=103) & ! eff.dat
            it_made, &
            int_abs2_a_plus_at_zl/lx, &
            loss_on_the_way_plus/lx, &
            int_abs2_a_plus_at_z0/lx, &
            sum_eff, &
            int_obrezki_z0, &
            int_obrezki_zl, &
            eff_tmp, &
            eff_tmp/(4.0d0*c3*sum_eff), &
            cdabs(a1(ix_out))

        eff_tmp_b = (loss_on_the_way_minus + int_abs2_a_minus_at_z0 - int_abs2_a_minus_at_zl_on_mir)/lx

        write (35, 108, err=103) & ! eff_b.dat
            it_made, &
            loss_on_the_way_minus/lx, &
            int_abs2_a_minus_at_zl_on_mir/lx, &
            int_abs2_a_minus_at_z0/lx, &
            eff_tmp_b

108     format(i, 4e17.8)

        write (33, 104, err=103) & ! eff_k.dat
            it_made, &
            int_abs2_a_plus_at_zl_k, &
            loss_on_the_way_plus_k, &
            int_abs2_a_plus_at_z0_k, &
            sum_eff, &
            int_obrezki_z0, &
            int_obrezki_zl, &
            eff_tmp_k, &
            eff_tmp_k/(4.0d0*c3*sum_eff), &
            cdabs(a1(ix_out))

104     format(i, 9e17.8)

        eff_tmp_k_b = loss_on_the_way_minus_k + int_abs2_a_minus_at_z0_k - int_abs2_a_minus_at_zl_on_mir_k

        write (335, 108, err=103) & ! eff_k_b.dat
            it_made, &
            loss_on_the_way_minus_k, &
            int_abs2_a_minus_at_zl_on_mir_k, &
            int_abs2_a_minus_at_z0_k, &
            eff_tmp_k_b

        write (53, 107, err=103) it_made, & ! eff_new.dat
            int_abs2_a_minus_at_z0_out_mir/lx, &
            int_abs2_a_plus_at_zl_out_mir/lx, &
            (loss_on_the_way_minus_k + loss_on_the_way_plus/lx), &
            obrezki, &
            4.0d0*c3*sum_eff, &
            loss_on_the_way_minus_k + (int_abs2_a_plus_at_zl_out_mir + &
                                       loss_on_the_way_plus + &
                                       int_abs2_a_minus_at_z0_out_mir)/lx + obrezki - &
            4.0d0*c3*sum_eff, &
            (loss_on_the_way_minus_k + (int_abs2_a_plus_at_zl_out_mir + &
                                        loss_on_the_way_plus + &
                                        int_abs2_a_minus_at_z0_out_mir)/lx + obrezki - &
             4.0d0*c3*sum_eff)/4.0d0/c3/sum_eff

        write (533, 107, err=103) it_made, & ! eff_new_k.dat
            int_abs2_a_minus_at_z0_out_mir_k, &
            int_abs2_a_plus_at_zl_out_mir_k, &
            (loss_on_the_way_minus_k + loss_on_the_way_plus_k), &
            obrezki, &
            4.0d0*c3*sum_eff, &
            loss_on_the_way_minus_k + (int_abs2_a_plus_at_zl_out_mir_k + &
                                       loss_on_the_way_plus_k + &
                                       int_abs2_a_minus_at_z0_out_mir_k + obrezki) - &
            4.0d0*c3*sum_eff, &
            (loss_on_the_way_minus_k + (int_abs2_a_plus_at_zl_out_mir_k + &
                                        loss_on_the_way_plus_k + &
                                        int_abs2_a_minus_at_z0_out_mir_k + obrezki) - &
             4.0d0*c3*sum_eff)/4.0d0/c3/sum_eff

107     format(i, 7e17.8)

        write (777, '(8e17.8)', err=103) & ! for_graphics.dat
            lz, &
            int_abs2_a_minus_at_z0_out_mir_k, &
            int_abs2_a_plus_at_zl_out_mir_k, &
            (loss_on_the_way_minus_k + loss_on_the_way_plus_k), &
            obrezki, &
            4.0d0*c3*sum_eff, &
            (int_abs2_a_minus_at_z0_out_mir_k + int_abs2_a_plus_at_zl_out_mir_k)/(4.0d0*c3*sum_eff), &
            delta

        !ix_out = 65
        if (recount == .false.) then
            write (*, '(a,\)') char(13)
                write(*,'(a,i6,a,f7.3,a,e17.8,a,e17.8,a)') 'iteration #', it_made, ': a(', x_out, ') = ', cdabs(a1(ix_out)), ' eff = ', sum_eff, ' carryover'
        else
            write (*, '(a,\)') char(13)
                write(*,'(a,i6,a,f7.3,a,e17.8,a,e17.8,a)') 'iteration #', it_made, ': a(', x_out, ') = ', cdabs(a1(ix_out)), ' eff = ', sum_eff, ' recount'
        end if

        !call fn_for_fortran_to_call(ptr)
    end do do_iter

    write (*, '(/,/)')

    !closing write files a(z=0,x) and a(z=l,x)
    close (1)
    close (3)
    close (53)
    close (221)
    close (222)
    close (788)
    if (amp_only == .false.) then
        close (2)
        close (8)
    end if
    if (thout == .true.) close (10)

    call cpu_time(finish_time)
    print *, 'computation time = ', finish_time - start_time, 'seconds'

    call write_result()
    call finish()

    print *, 'calculation finished.'
    pause
    stop
101 stop 'error of file open.'
102 stop 'error of file reading.'
103 stop 'error of file writing.'

end program elektron2dsi
