!    The file contains the following subrouines:

!          abq_UEL_3D_integrationpoints           - defines integration ponits for 3D continuum elements
!          abq_UEL_2D_integrationpoints           - defines integration points for 2D continuum elements
!          abq_UEL_1D_integrationpoints(n_points, n_nodes, xi, w)  = defines integration points for 1D line integral

      subroutine abq_UEL_3D_integrationpoints(n_points, n_nodes, xi, w)

      implicit none
      integer, intent(in) :: n_points
      integer, intent(in) :: n_nodes

      double precision, intent(out) :: xi(3,*)
      double precision, intent(out) :: w(*)

      integer :: i,j,k,n

      double precision x1D(4), w1D(4)

    !         Defines integration points and weights for 3D continuum elements

      if (n_nodes  == 4.or.n_nodes ==10 ) then   ! Tetrahedral elements
        if (n_points == 1) then
            xi(1,1) = 0.25D0
            xi(2,1) = 0.25D0
            xi(3,1) = 0.25D0
            w(1) = 1.D0/6.D0
        else if (n_points == 4) then
            xi(1,1) = 0.58541020
            xi(2,1) = 0.13819660
            xi(3,1) = xi(2,1)
            xi(1,2) = xi(2,1)
            xi(2,2) = xi(1,1)
            xi(3,2) = xi(2,1)
            xi(1,3) = xi(2,1)
            xi(2,3) = xi(2,1)
            xi(3,3) = xi(1,1)
            xi(1,4) = xi(2,1)
            xi(2,4) = xi(2,1)
            xi(3,4) = xi(2,1)
            w(1:4) = 1.D0/24.D0
        else if (n_points == 5) then
            xi(1,1) = 0.25d0
            xi(2,1) = 0.25d0
            xi(3,1) = 0.25d0
            xi(1,2) = 0.5d0
            xi(2,2) = 1.d0/6.d0
            xi(3,2) = 1.d0/6.d0
            xi(1,3) = 1.d0/6.d0
            xi(2,3) = 0.5d0
            xi(3,3) = 1.d0/6.d0
            xi(1,4) = 1.d0/6.d0
            xi(2,4) = 1.d0/6.d0
            xi(3,4) = 0.5d0
            xi(1,5) = 1.d0/6.d0
            xi(2,5) = 1.d0/6.d0
            xi(3,5) = 1.d0/6.d0
            w(1) = -4.d0/30.d0
            w(2:5) = 3.d0/40.d0
        else
            write(6,*) 'Incorrect # of int pts for tetrahedral element '
            write(6, *) ' called with ',n_points
            stop
        endif
      else if ( n_nodes == 8 .or. n_nodes == 20 ) then   ! 8 or 20 noded hexahedral elements
        if (n_points == 1) then
            xi(1,1) = 0.D0
            xi(2,1) = 0.D0
            xi(3,1) = 0.D0
            w(1) = 8.D0
        else if (n_points == 8) then
            x1D(1) = -0.5773502692
            x1D(2) =  0.5773502692
            do k = 1,2
                do j = 1,2
                    do i = 1,2
                        n = 4*(k-1) + 2*(j-1) + i
                        xi(1,n) = x1D(i)
                        xi(2,n) = x1D(j)
                        xi(3,n) = x1D(k)
                    end do
                end do
            end do
            w(1:8) = 1.D0
        else if (n_points == 27) then
            x1D(1) = -0.7745966692
            x1D(2) = 0.
            x1D(3) = 0.7745966692
            w1D(1) = 0.5555555555D0
            w1D(2) = 0.888888888D0
            w1D(3) = 0.55555555555D0
            do k = 1,3
                do j = 1,3
                    do i = 1,3
                        n = 9*(k-1) + 3*(j-1) + i
                        xi(1,n) = x1D(i)
                        xi(2,n) = x1D(j)
                        xi(3,n) = x1D(k)
                        w(n) = w1D(i)*w1D(j)*w1D(k)
                    end do
                end do
            end do
        else if (n_points == 64) then
            x1D(1) = .8611363115940526D+00
            x1D(2) = .3399810435848563D+00
            x1D(3) = -.3399810435848563D+00
            x1D(4) = -.8611363115940526D+00
            w1D(1) = .3478548451374538D+00
            w1D(2) = .6521451548625461D+00
            w1D(3) = .6521451548625461D+00
            w1D(4) = .3478548451374538D+00
            do k = 1,4
                do j = 1,4
                    do i = 1,4
                        n = 16*(k-1) + 4*(j-1) + i
                        xi(1,n) = x1D(i)
                        xi(2,n) = x1D(j)
                        xi(3,n) = x1D(k)
                        w(n) = w1D(i)*w1D(j)*w1D(k)
                    end do
                end do
            end do
        endif
      endif

      return

      end subroutine abq_UEL_3D_integrationpoints


      subroutine abq_UEL_2D_integrationpoints(n_points, n_nodes, xi, w)
      implicit none
      integer, intent(in) :: n_points
      integer, intent(in) :: n_nodes

      double precision, intent(out) :: xi(2,*)
      double precision, intent(out) :: w(*)

      integer :: i,j,k,n

      double precision :: cn,w1,w2,w11,w12,w22

    !         Defines integration points and weights for 2D continuum elements

      if ( n_points==1 ) then
        if ( n_nodes==4 .or. n_nodes==9 ) then    !     ---   4 or 9 noded quad
            xi(1, 1) = 0.D0
            xi(2, 1) = 0.D0
            w(1) = 4.D0
        else if ( n_nodes==3 .or. n_nodes==6 ) then !     ---   3 or 6 noded triangle
            xi(1, 1) = 1.D0/3.D0
            xi(2, 1) = 1.D0/3.D0
            w(1) = 1.D0/2.D0
        end if
      else if ( n_points==3 ) then
        xi(1, 1) = 0.5D0
        xi(2, 1) = 0.5D0
        w(1) = 1.D0/6.D0
        xi(1, 2) = 0.D0
        xi(2, 2) = 0.5D0
        w(2) = w(1)
        xi(1, 3) = 0.5D0
        xi(2, 3) = 0.D0
        w(3) = w(1)
      else if ( n_points==4 ) then
        if ( n_nodes==4 .or. n_nodes==8 .or. n_nodes==9 ) then
            !     2X2 GAUSS INTEGRATION POINTS FOR QUADRILATERAL
            !     43
            !     12
            cn = 0.5773502691896260D0
            xi(1, 1) = -cn
            xi(1, 2) = cn
            xi(1, 3) = cn
            xi(1, 4) = -cn
            xi(2, 1) = -cn
            xi(2, 2) = -cn
            xi(2, 3) = cn
            xi(2, 4) = cn
            w(1) = 1.D0
            w(2) = 1.D0
            w(3) = 1.D0
            w(4) = 1.D0
        else if ( n_nodes==3 .or. n_nodes==6 ) then
            !     xi integration points for triangle
            xi(1, 1) = 1.D0/3.D0
            xi(2, 1) = xi(1, 1)
            w(1) = -27.D0/96.D0
            xi(1, 2) = 0.6D0
            xi(2, 2) = 0.2D0
            w(2) = 25.D0/96.D0
            xi(1, 3) = 0.2D0
            xi(2, 3) = 0.6D0
            w(3) = w(2)
            xi(1, 4) = 0.2D0
            xi(2, 4) = 0.2D0
            w(4) = w(2)
        end if

      else if ( n_points==7 ) then
        ! Quintic integration for triangle
        xi(1,1) = 1.d0/3.d0
        xi(2,1) = xi(1,1)
        w(1) = 0.1125d0
        xi(1,2) = 0.0597158717d0
        xi(2,2) = 0.4701420641d0
        w(2) = 0.0661970763d0
        xi(1,3) = xi(2,2)
        xi(2,3) = xi(1,2)
        w(3) = w(2)
        xi(1,4) = xi(2,2)
        xi(2,4) = xi(2,2)
        w(4) = w(2)
        xi(1,5) = 0.7974269853d0
        xi(2,5) = 0.1012865073d0
        w(5) = 0.0629695902d0
        xi(1,6) = xi(2,5)
        xi(2,6) = xi(1,5)
        w(6) = w(5)
        xi(1,7) = xi(2,5)
        xi(2,7) = xi(2,5)
        w(7) = w(5)
      else if ( n_points==9 ) then
        !     3X3 GAUSS INTEGRATION POINTS
        !     789
        !     456
        !     123
        cn = 0.7745966692414830D0
        xi(1, 1) = -cn
        xi(1, 2) = 0.D0
        xi(1, 3) = cn
        xi(1, 4) = -cn
        xi(1, 5) = 0.D0
        xi(1, 6) = cn
        xi(1, 7) = -cn
        xi(1, 8) = 0.D0
        xi(1, 9) = cn
        xi(2, 1) = -cn
        xi(2, 2) = -cn
        xi(2, 3) = -cn
        xi(2, 4) = 0.D0
        xi(2, 5) = 0.D0
        xi(2, 6) = 0.D0
        xi(2, 7) = cn
        xi(2, 8) = cn
        xi(2, 9) = cn
        w1 = 0.5555555555555560D0
        w2 = 0.8888888888888890D0
        w11 = w1*w1
        w12 = w1*w2
        w22 = w2*w2
        w(1) = w11
        w(2) = w12
        w(3) = w11
        w(4) = w12
        w(5) = w22
        w(6) = w12
        w(7) = w11
        w(8) = w12
        w(9) = w11
      end if

      return

      end subroutine abq_UEL_2D_integrationpoints

      subroutine abq_UEL_1D_integrationpoints(n_points, n_nodes, xi, w)
      implicit none
      integer, intent(in) :: n_points
      integer, intent(in) :: n_nodes

      double precision, intent(out) :: xi(*)
      double precision, intent(out) :: w(*)

      integer :: i,j,k,n

      double precision x1D(4), w1D(4)



      select case ( n_points )
        case (2)
            xi(1) = .5773502691896257D+00
            xi(2) = -.5773502691896257D+00
            w(1) = .1000000000000000D+01
            w(2) = .1000000000000000D+01
            return
        case (3)
            xi(1) = 0.7745966692414834D+00
            xi(2) = .0000000000000000D+00
            xi(3) = -.7745966692414834D+00
            w(1) = .5555555555555556D+00
            w(2) = .8888888888888888D+00
            w(3) = .5555555555555556D+00
            return
        case (4)
            xi(1) = .8611363115940526D+00
            xi(2) = .3399810435848563D+00
            xi(3) = -.3399810435848563D+00
            xi(4) = -.8611363115940526D+00
            w(1) = .3478548451374538D+00
            w(2) = .6521451548625461D+00
            w(3) = .6521451548625461D+00
            w(4) = .3478548451374538D+00
            return
        case (5)
            xi(1) = .9061798459386640D+00
            xi(2) = .5384693101056831D+00
            xi(3) = .0000000000000000D+00
            xi(4) = -.5384693101056831D+00
            xi(5) = -.9061798459386640D+00
            w(1) = .2369268850561891D+00
            w(2) = .4786286704993665D+00
            w(3) = .5688888888888889D+00
            w(4) = .4786286704993665D+00
            w(5) = .2369268850561891D+00
            return
        case (6)
            xi(1) = .9324695142031521D+00
            xi(2) = .6612093864662645D+00
            xi(3) = .2386191860831969D+00
            xi(4) = -.2386191860831969D+00
            xi(5) = -.6612093864662645D+00
            xi(6) = -.9324695142031521D+00
            w(1) = .1713244923791703D+00
            w(2) = .3607615730481386D+00
            w(3) = .4679139345726910D+00
            w(4) = .4679139345726910D+00
            w(5) = .3607615730481386D+00
            w(6) = .1713244923791703D+00
            return
        case DEFAULT
            write(6,*)'Error in subroutine abq_UEL_1D_integrationpoints'
            write(6,*) ' Invalid number of integration points for 1D'
            write(6,*) ' n_points must be between 1 and 6'
            stop
      end select

      end subroutine ABQ_UEL_1D_integrationpoints
      