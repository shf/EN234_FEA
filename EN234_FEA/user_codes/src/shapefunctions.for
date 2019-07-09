!    The file contains the following subrouines:
!          abq_UEL_2D_shapefunctions              - defines shape functions for 2D continuum elements
!          abq_UEL_3D_shapefunctions              - defines shape functions for 3D continuum elements
!          abq_UEL_1D_shapefunctions              - defines shape functions for 1D continuum elements

      subroutine abq_UEL_2D_shapefunctions(xi,n_nodes,f,df)

      implicit none
      integer, intent(in) :: n_nodes
  
      double precision, intent(in) :: xi(2)
      double precision, intent(out) :: f(*)
      double precision, intent(out) :: df(9,2)
      double precision g1, g2, g3, dg1, dg2, dg3
      double precision h1, h2, h3, dh1, dh2, dh3
      double precision z,dzdp, dzdq

      if ( n_nodes==3 ) then        !     SHAPE FUNCTIONS FOR 3 NODED TRIANGLE
          f(1) = xi(1)
          f(2) = xi(2)
          f(3) = 1.D0 - xi(1) - xi(2)
          df(1, 1) = 1.D0
          df(1, 2) = 0.D0
          df(2, 1) = 0.D0
          df(2, 2) = 1.D0
          df(3, 1) = -1.D0
          df(3, 2) = -1.D0
      else if ( n_nodes==4 ) then
          !     SHAPE FUNCTIONS FOR 4 NODED QUADRILATERAL
          !     43
          !     12
          g1 = 0.5D0*(1.D0 - xi(1))
          g2 = 0.5D0*(1.D0 + xi(1))
          h1 = 0.5D0*(1.D0 - xi(2))
          h2 = 0.5D0*(1.D0 + xi(2))
          f(1) = g1*h1
          f(2) = g2*h1
          f(3) = g2*h2
          f(4) = g1*h2
          dg1 = -0.5D0
          dg2 = 0.5D0
          dh1 = -0.5D0
          dh2 = 0.5D0
          df(1, 1) = dg1*h1
          df(2, 1) = dg2*h1
          df(3, 1) = dg2*h2
          df(4, 1) = dg1*h2
          df(1, 2) = g1*dh1
          df(2, 2) = g2*dh1
          df(3, 2) = g2*dh2
          df(4, 2) = g1*dh2
      else if ( n_nodes==6 ) then
          !     SHAPE FUNCTIONS FOR 6 NODED TRIANGLE
          !          3
          !       6      5
          !     1    4     2
          !     P = L1
          !     Q = L2
          !     Z = 1 - P - Q = L3
          z = 1.D0 - xi(1) - xi(2)
          f(1) = (2.D0*xi(1) - 1.D0)*xi(1)
          f(2) = (2.D0*xi(2) - 1.D0)*xi(2)
          f(3) = (2.D0*z - 1.D0)*z
          f(4) = 4.D0*xi(1)*xi(2)
          f(5) = 4.D0*xi(2)*z
          f(6) = 4.D0*xi(1)*z
          dzdp = -1.D0
          dzdq = -1.D0
          df(1, 1) = 4.D0*xi(1) - 1.D0
          df(2, 1) = 0.D0
          df(3, 1) = 4.D0*z*dzdp - dzdp
          df(4, 1) = 4.D0*xi(2)
          df(5, 1) = 4.D0*xi(2)*dzdp
          df(6, 1) = 4.D0*z + 4.D0*xi(1)*dzdp
          df(1, 2) = 0.D0
          df(2, 2) = 4.D0*xi(2) - 1.D0
          df(3, 2) = 4.D0*z*dzdq - dzdq
          df(4, 2) = 4.D0*xi(1)
          df(5, 2) = 4.D0*z + 4.D0*xi(2)*dzdq
          df(6, 2) = 4.D0*xi(1)*dzdq
      else if ( n_nodes==8 ) then
          !     SHAPE FUNCTIONS FOR 8 NODED SERENDIPITY ELEMENT
           f(1) = -0.25*(1.-xi(1))*(1.-xi(2))*(1.+xi(1)+xi(2));
           f(2) = 0.25*(1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)-1.);
           f(3) = 0.25*(1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)-1.);
           f(4) = 0.25*(1.-xi(1))*(1.+xi(2))*(xi(2)-xi(1)-1.);
           f(5) = 0.5*(1.-xi(1)*xi(1))*(1.-xi(2));
           f(6) = 0.5*(1.+xi(1))*(1.-xi(2)*xi(2));
           f(7) = 0.5*(1.-xi(1)*xi(1))*(1.+xi(2));
           f(8) = 0.5*(1.-xi(1))*(1.-xi(2)*xi(2));
           df(1,1) = 0.25*(1.-xi(2))*(2.*xi(1)+xi(2));
           df(1,2) = 0.25*(1.-xi(1))*(xi(1)+2.*xi(2));
           df(2,1) = 0.25*(1.-xi(2))*(2.*xi(1)-xi(2));
           df(2,2) = 0.25*(1.+xi(1))*(2.*xi(2)-xi(1));
           df(3,1) = 0.25*(1.+xi(2))*(2.*xi(1)+xi(2));
           df(3,2) = 0.25*(1.+xi(1))*(2.*xi(2)+xi(1));
           df(4,1) = 0.25*(1.+xi(2))*(2.*xi(1)-xi(2));
           df(4,2) = 0.25*(1.-xi(1))*(2.*xi(2)-xi(1));
           df(5,1) = -xi(1)*(1.-xi(2));
           df(5,2) = -0.5*(1.-xi(1)*xi(1));
           df(6,1) = 0.5*(1.-xi(2)*xi(2));
           df(6,2) = -(1.+xi(1))*xi(2);
           df(7,1) = -xi(1)*(1.+xi(2));
           df(7,2) = 0.5*(1.-xi(1)*xi(1));
           df(8,1) = -0.5*(1.-xi(2)*xi(2));
           df(8,2) = -(1.-xi(1))*xi(2);
      else if ( n_nodes==9 ) then
          !     SHAPE FUNCTIONS FOR 9 NODED LAGRANGIAN ELEMENT
          !     789
          !     456
          !     123
          g1 = -.5D0*xi(1)*(1.D0 - xi(1))
          g2 = (1.D0 - xi(1))*(1.D0 + xi(1))
          g3 = .5D0*xi(1)*(1.D0 + xi(1))
          h1 = -.5D0*xi(2)*(1.D0 - xi(2))
          h2 = (1.D0 - xi(2))*(1.D0 + xi(2))
          h3 = .5D0*xi(2)*(1.D0 + xi(2))
          dg1 = xi(1) - 0.5d0
          dg2 = -2.d0*xi(1)
          dg3 = xi(1) + 0.5d0
          dh1 = xi(2)-0.5d0
          dh2 = -2.d0*xi(2)
          dh3 = xi(2) + 0.5d0
          f(1) = g1*h1
          f(2) = g2*h1
          f(3) = g3*h1
          f(4) = g1*h2
          f(5) = g2*h2
          f(6) = g3*h2
          f(7) = g1*h3
          f(8) = g2*h3
          f(9) = g3*h3
          df(1,1) = dg1*h1
          df(1,2) = g1*dh1
          df(2,1) = dg2*h1
          df(2,2) = g2*dh1
          df(3,1) = dg3*h1
          df(3,2) = g3*dh1
          df(4,1) = dg1*h2
          df(4,2) = g1*dh2
          df(5,1) = dg2*h2
          df(5,2) = g2*dh2
          df(6,1) = dg3*h2
          df(6,2) = g3*dh2
          df(7,1) = dg1*h3
          df(7,2) = g1*dh3
          df(8,1) = dg2*h3
          df(8,2) = g2*dh3
          df(9,1) = dg3*h3
          df(9,2) = g3*dh3
      end if

      end subroutine abq_UEL_2D_shapefunctions

      subroutine abq_UEL_1D_shapefunctions(xi,n_nodes,f,df)
	  
        implicit none
        integer, intent(in) :: n_nodes
      
        double precision, intent(in) :: xi
        double precision, intent(out) :: f(*)
        double precision, intent(out) :: df(3)
      
            if (n_nodes==2) then
                f(1) = 0.5d0*(1.d0-xi)
                f(2) = 0.5d0*(1.d0+xi)
                df(1) = -0.5d0
                df(2) =  0.5d0
            else if (n_nodes==3) then
                f(1) = -0.5*xi*(1.-xi)
                f(2) =  0.5*xi*(1.+xi)
                f(3) = (1.-xi)*(1.+xi)
                df(1) = -0.5+xi
                df(2) =  0.5+xi
                df(3) = -2.d0*xi
            endif

        end subroutine abq_UEL_1D_shapefunctions

      subroutine abq_UEL_3D_shapefunctions(xi,n_nodes,f,df)
          implicit none
          integer, intent(in) :: n_nodes
    
          double precision, intent(in) :: xi(3)
          double precision, intent(out) :: f(20)
          double precision, intent(out) :: df(20,3)
          double precision xi4

!   Defines shape functions for 3D continuum elements

      if (n_nodes == 4) then
        f(1) = xi(1)
        f(2) = xi(2)
        f(3) = xi(3)
        f(4) = 1.-xi(1)-xi(2)-xi(3)
        df(1,1) = 1.
        df(2,2) = 1.
        df(3,3) = 1.
        df(4,1) = -1.
        df(4,2) = -1.
        df(4,3) = -1.
      else if (n_nodes == 10) then
        xi4 = 1.D0-xi(1)-xi(2)-xi(3)
        f(1) = (2.*xi(1)-1.)*xi(1)
        f(2) = (2.*xi(2)-1.)*xi(2)
        f(3) = (2.*xi(3)-1.)*xi(3)
        f(4) = (2.*xi4-1.)*xi4
        f(5) = 4.*xi(1)*xi(2)
        f(6) = 4.*xi(2)*xi(3)
        f(7) = 4.*xi(3)*xi(1)
        f(8) = 4.*xi(1)*xi4
        f(9) = 4.*xi(2)*xi4
        f(10) = 4.*xi(3)*xi4
        df(1,1) = (4.*xi(1)-1.)
        df(2,2) = (4.*xi(2)-1.)
        df(3,3) = (4.*xi(3)-1.)
        df(4,1) = -(4.*xi4-1.)
        df(4,2) = -(4.*xi4-1.)
        df(4,3) = -(4.*xi4-1.)
        df(5,1) = 4.*xi(2)
        df(5,2) = 4.*xi(1)
        df(6,2) = 4.*xi(3)
        df(6,3) = 4.*xi(2)
        df(7,1) = 4.*xi(3)
        df(7,3) = 4.*xi(1)
        df(8,1) = 4.*(xi4-xi(1))
        df(8,2) = -4.*xi(1)
        df(8,3) = -4.*xi(1)
        df(9,1) = -4.*xi(2)
        df(9,2) = 4.*(xi4-xi(2))
        df(9,3) = -4.*xi(2)
        df(10,1) = -4.*xi(3)*xi4
        df(10,2) = -4.*xi(3)
        df(10,3) = 4.*(xi4-xi(3))
      else if (n_nodes == 8) then
        f(1) = (1.-xi(1))*(1.-xi(2))*(1.-xi(3))/8.
        f(2) = (1.+xi(1))*(1.-xi(2))*(1.-xi(3))/8.
        f(3) = (1.+xi(1))*(1.+xi(2))*(1.-xi(3))/8.
        f(4) = (1.-xi(1))*(1.+xi(2))*(1.-xi(3))/8.
        f(5) = (1.-xi(1))*(1.-xi(2))*(1.+xi(3))/8.
        f(6) = (1.+xi(1))*(1.-xi(2))*(1.+xi(3))/8.
        f(7) = (1.+xi(1))*(1.+xi(2))*(1.+xi(3))/8.
        f(8) = (1.-xi(1))*(1.+xi(2))*(1.+xi(3))/8.
        df(1,1) = -(1.-xi(2))*(1.-xi(3))/8.
        df(1,2) = -(1.-xi(1))*(1.-xi(3))/8.
        df(1,3) = -(1.-xi(1))*(1.-xi(2))/8.
        df(2,1) = (1.-xi(2))*(1.-xi(3))/8.
        df(2,2) = -(1.+xi(1))*(1.-xi(3))/8.
        df(2,3) = -(1.+xi(1))*(1.-xi(2))/8.
        df(3,1) = (1.+xi(2))*(1.-xi(3))/8.
        df(3,2) = (1.+xi(1))*(1.-xi(3))/8.
        df(3,3) = -(1.+xi(1))*(1.+xi(2))/8.
        df(4,1) = -(1.+xi(2))*(1.-xi(3))/8.
        df(4,2) = (1.-xi(1))*(1.-xi(3))/8.
        df(4,3) = -(1.-xi(1))*(1.+xi(2))/8.
        df(5,1) = -(1.-xi(2))*(1.+xi(3))/8.
        df(5,2) = -(1.-xi(1))*(1.+xi(3))/8.
        df(5,3) = (1.-xi(1))*(1.-xi(2))/8.
        df(6,1) = (1.-xi(2))*(1.+xi(3))/8.
        df(6,2) = -(1.+xi(1))*(1.+xi(3))/8.
        df(6,3) = (1.+xi(1))*(1.-xi(2))/8.
        df(7,1) = (1.+xi(2))*(1.+xi(3))/8.
        df(7,2) = (1.+xi(1))*(1.+xi(3))/8.
        df(7,3) = (1.+xi(1))*(1.+xi(2))/8.
        df(8,1) = -(1.+xi(2))*(1.+xi(3))/8.
        df(8,2) = (1.-xi(1))*(1.+xi(3))/8.
        df(8,3) = (1.-xi(1))*(1.+xi(2))/8.
      else if (n_nodes == 20) then
        f(1)=(1.-xi(1))*(1.-xi(2))*(1.-xi(3))*(-xi(1)-xi(2)-xi(3)-2.)/8.
        f(2)=(1.+xi(1))*(1.-xi(2))*(1.-xi(3))*(xi(1)-xi(2)-xi(3)-2.)/8.
        f(3)=(1.+xi(1))*(1.+xi(2))*(1.-xi(3))*(xi(1)+xi(2)-xi(3)-2.)/8.
        f(4)=(1.-xi(1))*(1.+xi(2))*(1.-xi(3))*(-xi(1)+xi(2)-xi(3)-2.)/8.
        f(5)=(1.-xi(1))*(1.-xi(2))*(1.+xi(3))*(-xi(1)-xi(2)+xi(3)-2.)/8.
        f(6)=(1.+xi(1))*(1.-xi(2))*(1.+xi(3))*(xi(1)-xi(2)+xi(3)-2.)/8.
        f(7)=(1.+xi(1))*(1.+xi(2))*(1.+xi(3))*(xi(1)+xi(2)+xi(3)-2.)/8.
        f(8)=(1.-xi(1))*(1.+xi(2))*(1.+xi(3))*(-xi(1)+xi(2)+xi(3)-2.)/8.
        f(9) = (1.-xi(1)**2.)*(1.-xi(2))*(1.-xi(3))/4.
        f(10) = (1.+xi(1))*(1.-xi(2)**2.)*(1.-xi(3))/4.
        f(11) = (1.-xi(1)**2.)*(1.+xi(2))*(1.-xi(3))/4.
        f(12) = (1.-xi(1))*(1.-xi(2)**2.)*(1.-xi(3))/4.
        f(13) = (1.-xi(1)**2.)*(1.-xi(2))*(1.+xi(3))/4.
        f(14) = (1.+xi(1))*(1.-xi(2)**2.)*(1.+xi(3))/4.
        f(15) = (1.-xi(1)**2.)*(1.+xi(2))*(1.+xi(3))/4.
        f(16) = (1.-xi(1))*(1.-xi(2)**2.)*(1.+xi(3))/4.
        f(17) = (1.-xi(1))*(1.-xi(2))*(1.-xi(3)**2.)/4.
        f(18) = (1.+xi(1))*(1.-xi(2))*(1.-xi(3)**2.)/4.
        f(19) = (1.+xi(1))*(1.+xi(2))*(1.-xi(3)**2.)/4.
        f(20) = (1.-xi(1))*(1.+xi(2))*(1.-xi(3)**2.)/4.
        df(1,1) = (-(1.-xi(2))*(1.-xi(3))*(-xi(1)-xi(2)-xi(3)-2.)
     1           -(1.-xi(1))*(1.-xi(2))*(1.-xi(3)))/8.
        df(1,2) = (-(1.-xi(1))*(1.-xi(3))*(-xi(1)-xi(2)-xi(3)-2.)
     1           -(1.-xi(1))*(1.-xi(2))*(1.-xi(3)))/8.
        df(1,3) = (-(1.-xi(1))*(1.-xi(2))*(-xi(1)-xi(2)-xi(3)-2.)
     1           -(1.-xi(1))*(1.-xi(2))*(1.-xi(3)))/8.

        df(2,1) = ((1.-xi(2))*(1.-xi(3))*(xi(1)-xi(2)-xi(3)-2.)
     1           +(1.+xi(1))*(1.-xi(2))*(1.-xi(3)))/8.
        df(2,2) = (-(1.+xi(1))*(1.-xi(3))*(xi(1)-xi(2)-xi(3)-2.)
     1          -(1.+xi(1))*(1.-xi(2))*(1.-xi(3)))/8.
        df(2,3) = (-(1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)-xi(3)-2.)
     1           -(1.+xi(1))*(1.-xi(2))*(1.-xi(3)))/8.

        df(3,1) = ((1.+xi(2))*(1.-xi(3))*(xi(1)+xi(2)-xi(3)-2.)
     1           +(1.+xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
        df(3,2) = ((1.+xi(1))*(1.-xi(3))*(xi(1)+xi(2)-xi(3)-2.)
     1           +(1.+xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
        df(3,3) = (-(1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)-xi(3)-2.)
     1           -(1.+xi(1))*(1.+xi(2))*(1.-xi(3)))/8.

        df(4,1) = (-(1.+xi(2))*(1.-xi(3))*(-xi(1)+xi(2)-xi(3)-2.)
     1           -(1.-xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
        df(4,2) = ((1.-xi(1))*(1.-xi(3))*(-xi(1)+xi(2)-xi(3)-2.)
     1            +(1.-xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
        df(4,3) = (-(1.-xi(1))*(1.+xi(2))*(-xi(1)+xi(2)-xi(3)-2.)
     1           -(1.-xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
        df(5,1) = (-(1.-xi(2))*(1.+xi(3))*(-xi(1)-xi(2)+xi(3)-2.)
     1           -(1.-xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
        df(5,2) = (-(1.-xi(1))*(1.+xi(3))*(-xi(1)-xi(2)+xi(3)-2.)
     1           -(1.-xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
        df(5,3) = ((1.-xi(1))*(1.-xi(2))*(-xi(1)-xi(2)+xi(3)-2.)
     1           +(1.-xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
        df(6,1) = ((1.-xi(2))*(1.+xi(3))*(xi(1)-xi(2)+xi(3)-2.)
     1           +(1.+xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
        df(6,2) = (-(1.+xi(1))*(1.+xi(3))*(xi(1)-xi(2)+xi(3)-2.)
     1           -(1.+xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
        df(6,3) = ((1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)+xi(3)-2.)
     1           +(1.+xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
        df(7,1) = ((1.+xi(2))*(1.+xi(3))*(xi(1)+xi(2)+xi(3)-2.)
     1           +(1.+xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
        df(7,2) = ((1.+xi(1))*(1.+xi(3))*(xi(1)+xi(2)+xi(3)-2.)
     1           +(1.+xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
        df(7,3) = ((1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)+xi(3)-2.)
     1           +(1.+xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
        df(8,1) = (-(1.+xi(2))*(1.+xi(3))*(-xi(1)+xi(2)+xi(3)-2.)
     1           -(1.-xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
        df(8,2) = ((1.-xi(1))*(1.+xi(3))*(-xi(1)+xi(2)+xi(3)-2.)
     1           +(1.-xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
        df(8,3) = ((1.-xi(1))*(1.+xi(2))*(-xi(1)+xi(2)+xi(3)-2.)
     1           +(1.-xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
        df(9,1)  = -2.*xi(1)*(1.-xi(2))*(1.-xi(3))/4.
        df(9,2)  = -(1.-xi(1)**2.)*(1.-xi(3))/4.
        df(9,3)  = -(1.-xi(1)**2.)*(1.-xi(2))/4.
        df(10,1)  = (1.-xi(2)**2.)*(1.-xi(3))/4.
        df(10,2)  = -2.*xi(2)*(1.+xi(1))*(1.-xi(3))/4.
        df(10,3)  = -(1.-xi(2)**2.)*(1.+xi(1))/4.
        df(11,1)  = -2.*xi(1)*(1.-xi(2))*(1.-xi(3))/4.
        df(11,2)  = -(1.-xi(1)**2.)*(1.-xi(3))/4.
        df(11,3)  = -(1.-xi(1)**2.)*(1.-xi(2))/4.
        df(12,1)  = -(1.-xi(2)**2.)*(1.-xi(3))/4.
        df(12,2)  = -2.*xi(2)*(1.-xi(1))*(1.-xi(3))/4.
        df(12,3)  = -(1.-xi(2)**2.)*(1.-xi(1))/4.
        df(13,1)  = -2.*xi(1)*(1.-xi(2))*(1.+xi(3))/4.
        df(13,2)  = -(1.-xi(1)**2.)*(1.+xi(3))/4.
        df(13,3)  = (1.-xi(1)**2.)*(1.-xi(2))/4.
        df(14,1)  = (1.-xi(2)**2.)*(1.+xi(3))/4.
        df(14,2)  = -2.*xi(2)*(1.+xi(1))*(1.+xi(3))/4.
        df(14,3)  = (1.-xi(2)**2.)*(1.+xi(1))/4.
        df(15,1)  = 2.*xi(1)*(1.+xi(2))*(1.+xi(3))/4.
        df(15,2)  = (1.-xi(1)**2.)*(1.+xi(3))/4.
        df(15,3)  = (1.-xi(1)**2.)*(1.+xi(2))/4.
        df(16,1)  = -(1.-xi(2)**2.)*(1.+xi(3))/4.
        df(16,2)  = -2.*xi(2)*(1.-xi(1))*(1.+xi(3))/4.
        df(16,3)  = (1.-xi(2)**2.)*(1.-xi(1))/4.
        df(17,1) = -(1.-xi(2))*(1.-xi(3)**2.)/4.
        df(17,2) = -(1.-xi(1))*(1.-xi(3)**2.)/4.
        df(17,3) = -xi(3)*(1.-xi(1))*(1.-xi(2))/2.
        df(18,1) = (1.-xi(2))*(1.-xi(3)**2.)/4.
        df(18,2) = -(1.+xi(1))*(1.-xi(3)**2.)/4.
        df(18,3) = -xi(3)*(1.+xi(1))*(1.-xi(2))/2.
        df(19,1) = (1.+xi(2))*(1.-xi(3)**2.)/4.
        df(19,2) = (1.+xi(1))*(1.-xi(3)**2.)/4.
        df(19,3) = -xi(3)*(1.+xi(1))*(1.+xi(2))/2.
        df(20,1) = -(1.+xi(2))*(1.-xi(3)**2.)/4.
        df(20,2) = (1.-xi(1))*(1.-xi(3)**2.)/4.
        df(20,3) = -xi(3)*(1.-xi(1))*(1.+xi(2))/2.
      endif


      end subroutine abq_UEL_3D_shapefunctions
