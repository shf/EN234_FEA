!
!    ABAQUS format UEL subroutine
!
!    This file is compatible with both EN234_FEA and ABAQUS/Standard
!
!    The example implements a standard fully integrated 3D linear elastic continuum element
!
!=========================== ABAQUS format user element subroutine ===================

      SUBROUTINE UEL_PH(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3     LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
    !
      INCLUDE 'ABA_PARAM.INC'
    !
    !
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1   SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2   DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3   JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4   PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

    !
    !       Variables that must be computed in this routine
    !       RHS(i)                     Right hand side vector.  In EN234_FEA the dimensions are always RHS(MLVARX,1)
    !       AMATRX(i,j)                Stiffness matrix d RHS(i)/ d DU(j)
    !       SVARS(1:NSVARS)            Element state variables.  Must be updated in this routine
    !       ENERGY(1:8)
    !                                  Energy(1) Kinetic Energy
    !                                  Energy(2) Elastic Strain Energy
    !                                  Energy(3) Creep Dissipation
    !                                  Energy(4) Plastic Dissipation
    !                                  Energy(5) Viscous Dissipation
    !                                  Energy(6) Artificial strain energy
    !                                  Energy(7) Electrostatic energy
    !                                  Energy(8) Incremental work done by loads applied to the element
    !       PNEWDT                     Allows user to control ABAQUS time increments.
    !                                  If PNEWDT<1 then time step is abandoned and computation is restarted with
    !                                  a time increment equal to PNEWDT*DTIME
    !                                  If PNEWDT>1 ABAQUS may increase the time increment by a factor PNEWDT
    !
    !       Variables provided for information
    !       NDOFEL                     Total # DOF for the element
    !       NRHS                       Dimension variable
    !       NSVARS                     Total # element state variables
    !       PROPS(1:NPROPS)            User-specified properties of the element
    !       NPROPS                     No. properties
    !       JPROPS(1:NJPROPS)          Integer valued user specified properties for the element
    !       NJPROPS                    No. integer valued properties
    !       COORDS(i,N)                ith coordinate of Nth node on element
    !       MCRD                       Maximum of (# coords,minimum of (3,#DOF)) on any node
    !       U                          Vector of DOF at the end of the increment
    !       DU                         Vector of DOF increments
    !       V                          Vector of velocities (defined only for implicit dynamics)
    !       A                          Vector of accelerations (defined only for implicit dynamics)
    !       JTYPE                      Integer identifying element type (the number n in the Un specification in the input file)
    !       TIME(1:2)                  TIME(1)   Current value of step time
    !                                  TIME(2)   Total time
    !       DTIME                      Time increment
    !       KSTEP                      Current step number (always 1 in EN234_FEA)
    !       KINC                       Increment number
    !       JELEM                      User assigned element number in ABAQUS (internally assigned in EN234_FEA)
    !       PARAMS(1:3)                Time increment parameters alpha, beta, gamma for implicit dynamics
    !       NDLOAD                     Number of user-defined distributed loads defined for this element
    !       JDLTYP(1:NDLOAD)           Integers n defining distributed load types defined as Un or (if negative) UnNU in input file
    !       ADLMAG(1:NDLOAD)           Distributed load magnitudes
    !       DDLMAG(1:NDLOAD)           Increment in distributed load magnitudes
    !       PREDEF(1:2,1:NPREDF,1:NNODE)   Predefined fields.
    !       PREDEF(1,...)              Value of predefined field
    !       PREDEF(2,...)              Increment in predefined field
    !       PREDEF(1:2,1,k)            Value of temperature/temperature increment at kth node
    !       PREDEF(1:2,2:NPREDF,k)     Value of user defined field/field increment at kth node (not used in EN234FEA)
    !       NPREDF                     Number of predefined fields (1 for en234FEA)
    !       LFLAGS                     Control variable
    !       LFLAGS(1)                  Defines procedure type
    !       LFLAGS(2)                  0 => small displacement analysis  1 => Large displacement (NLGEOM option)
    !       LFLAGS(3)                   1 => Subroutine must return both RHS and AMATRX (always true in EN234FEA)
    !                                   2 => Subroutine must return stiffness AMATRX = -dF/du
    !                                   3 => Subroutine must return daming matrix AMATRX = -dF/dudot
    !                                   4 => Subroutine must return mass matrix AMATRX = -dF/duddot
    !                                   5 => Define the RHS only
    !                                   6 => Define the mass matrix for the initial acceleration calculation
    !                                   100 => Define perturbation quantities for output
    !       LFLAGS(4)                   0 => General step   1 => linear perturbation step
    !       LFLAGS(5)                   0 => current approximation to solution based on Newton correction; 1 => based on extrapolation
    !       MLVARX                      Dimension variable (equal to NDOFEL in EN234FEA)
    !       PERIOD                      Time period of the current step
    !
    !
    ! Local Variables
      integer      :: i,j,n_points,kint, nfacenodes, ipoin, ksize
      integer      :: face_node_list(3)                       ! List of nodes on an element face
    !
      double precision  ::  xi(2,9)                          ! Area integration points
      double precision  ::  xi_ph(2,4)                       ! phase field integration points
      double precision  ::  w(9)                             ! Area integration weights
      double precision  ::  w_ph(4)                          ! Phase field weights
      double precision  ::  N(9)                             ! 2D shape functions
      double precision  ::  N_ph(4)                          ! phase field shape functions
      double precision  ::  dNdxi(9,2)                       ! 2D shape function derivatives
      double precision  ::  dNdxi_ph(9,2)                    ! 2D shape function derivatives for phase field
      double precision  ::  dNdx(9,2)                        ! Spatial derivatives
      double precision  ::  dNdx_ph(9,2)                     ! Spatial derivatives for phase field
      double precision  ::  dxdxi(2,2)                       ! Derivative of spatial coords wrt normalized coords
      double precision  ::  dxdxi_ph(2,2)                    ! Derivative of spatial coords wrt normalized coords for phase field
      double precision  ::  dxidx(2,2), determinant, det0    ! Jacobian inverse and determinant
      double precision  ::  dxidx_ph(2,2), det_ph

    !   Variables below are for computing integrals over element faces
      double precision  ::  face_coords(2,3)                  ! Coords of nodes on an element face
      double precision  ::  xi1(6)                            ! 1D integration points
      double precision  ::  w1(6)                              ! Integration weights
      double precision  ::  N1(3)                             ! 1D shape functions
      double precision  ::  dNdxi1(3)                         ! 1D shape function derivatives
      double precision  ::  norm(2)                           ! Normal to an element face
      double precision  ::  dxdxi1(2)                         ! Derivative of 1D spatial coord wrt normalized areal coord
    !
      double precision  ::  strain(9)                         ! Strain vector contains [e11, e22, 2e12, mu, c, dmu/dx1, dmu/dx2, dc/dx1, dc/dx2]
      double precision  ::  dstrain(9)                        ! dStrain vector
      double precision  ::  q(9)                              ! node residual
      double precision  ::  elstrain(4)                        ! dStrain vector
      double precision  ::  elstress(4)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
      double precision  ::  D(9,9)                            ! D*(strain)  (NOTE FACTOR OF 2 in shear strain)
      double precision  ::  Del(4,4)                          ! elastic D
      double precision  ::  B(9,24)                           ! strain = B*(dof_total)
    !  double precision  ::  ktemp(22,22)                      ! Temporary stiffness (for incompatible mode elements)
    !  double precision  ::  rhs_temp(22)                      ! Temporary RHS vector (for incompatible mode elements)
    !  double precision  ::  kuu(18,18)                        ! Upper diagonal stiffness
    !  double precision  ::  kaa(4,4),kaainv(4,4)              ! Lower diagonal stiffness
    !  double precision  ::  kau(4,18)                         ! Lower quadrant of stiffness
    !  double precision  ::  kua(18,4)                         ! Upper quadrant of stiffness
    !  double precision  ::  alpha(4)                          ! Internal DOF for incompatible mode element
    !  double precision  ::  xm(2), wm(1)                         ! Integration point and weight for incompatible mode
      double precision  ::  E, xnu, Om, WC, kappa, dc, theta   ! Material properties
      double precision  ::  D44, D11, D12, con                 ! temporary
      double precision  ::  delta(4)                           ! delta

    !
    !     Example ABAQUS UEL implemeting 2D phase field model
    !     El props are:

    !     PROPS(1)         Young's modulus
    !     PROPS(2)         Poisson's ratio
    !     PROPS(3)         Volume fraction
    !     PROPS(4)         Energy of deviatoin
    !     PROPS(5)         Interface energy
    !     PROPS(6)         Diffusion coefficient
    !     PROPS(7)         Numerical scheme
       
    !  if (NNODE == 3) n_points = 1              ! Linear triangle
    !  if (NNODE == 4) n_points = 4               ! Linear rectangle
    !  if (NNODE == 6) n_points = 4              ! Quadratic triangle
      if (NNODE == 8) then
        n_points = 4               ! Serendipity rectangle
      else
        print*, 'Number of nodes is not correct'
      endif
    !  if (NNODE == 9) n_points = 9             ! Quadratic rect

      call abq_UEL_2D_integrationpoints(n_points, NNODE, xi, w)
      call abq_UEL_2D_integrationpoints(4, 4, xi_ph, w_ph)

	  
      if (MLVARX<2*NNODE) then
        write(6,*) ' Error in abaqus UEL '
        write(6,*) ' Variable MLVARX must exceed 2*NNODE'
        write(6,*) ' MLVARX = ',MLVARX,' NNODE = ',NNODE
        stop
      endif
	  
      RHS(1:MLVARX,1) = 0.d0
      AMATRX(1:NDOFEL,1:NDOFEL) = 0.d0

      SVARS(1:NSVARS) = 0.d0

      delta = (/1.d0,1.d0,1.d0,0.d0/)

    !  alpha = 0.d0
    !  kuu = 0.d0
    !  kau = 0.d0
    !  kaa = 0.d0
    !  kua = 0.d0
    !  rhs_temp = 0.D0
	  
      E = PROPS(1)
      xnu = PROPS(2)
      Om = PROPS(3)
      WC = PROPS(4)
      kappa = PROPS(5)
      dc = PROPS(6)
      theta = PROPS(7)

      h = 1.D0 ! thickness

      d11 = (1.D0-xnu)*E/( (1.D0+xnu)*(1.D0-2.D0*xnu) )
      d12 = xnu*E/( (1.D0+xnu)*(1.D0-2.D0*xnu) )
      d44 = (0.5D0 - xnu)*E/( (1.D0+xnu)*(1.D0-2.D0*xnu) )
      
      Del = 0.D0
      D = 0.D0
      
      Del(1:3,1:3) = d12
      Del(1,1) = d11
      Del(2,2) = d11
      Del(3,3) = d11
      Del(4,4) = d44

      D(1:2,1:2) = Del(1:2,1:2)
      D(3,3) = Del(4,4)


      D(4,1) = -Om*(Del(1, 1) + Del(1, 2) + Del(1, 3))/3.D0
      D(4,2) = -Om*(Del(2, 1) + Del(2, 2) + Del(2, 3))/3.D0
      D(1,5) = -Om*(Del(1, 1) + Del(1, 2) + Del(1, 3))/3.D0
      D(2,5) = -Om*(Del(2, 1) + Del(2, 2) + Del(2, 3))/3.D0
      D(4,4) = 1.D0
      D(5,5) = 1.D0 / DTIME
      D(6,8) = -kappa
      D(7,9) = -kappa
      D(8,6) = theta*dc
      D(9,7) = theta*dc
      
      ENERGY(1:8) = 0.d0
    
!      JTYPE = 2

      if (JTYPE==1) then

      do kint = 1, n_points
        strain = 0.d0
        dstrain = 0.d0
        elstress = 0.d0
        elstrain = 0.D0
        q = 0.D0
        con = 0.D0

        call abq_UEL_2D_shapefunctions(xi(1:2,kint),NNODE,N,dNdxi)
        dxdxi = matmul(coords(1:2,1:NNODE),dNdxi(1:NNODE,1:2))
        call abq_UEL_invert2d(dxdxi,dxidx,determinant)
        dNdx(1:NNODE,1:2) = matmul(dNdxi(1:NNODE,1:2),dxidx)

        call abq_UEL_2D_shapefunctions(xi_ph(1:2,kint),4,N_ph,dNdxi_ph)
        dxdxi_ph = matmul(coords(1:2,1:4),dNdxi(1:4,1:2))
        call abq_UEL_invert2d(dxdxi_ph,dxidx_ph,det_ph)
        dNdx_ph(1:4,1:2) = matmul(dNdxi(1:4,1:2),dxidx_ph) 

        B = 0.D0
        B(1,1)  = dNdx(1,1)
        B(1,5)  = dNdx(2,1)
        B(1,9)  = dNdx(3,1)
        B(1,13) = dNdx(4,1)
        B(1,17) = dNdx(5,1)
        B(1,19) = dNdx(6,1)
        B(1,21) = dNdx(7,1)
        B(1,23) = dNdx(8,1)

        B(2,2)  = dNdx(1,2)
        B(2,6)  = dNdx(2,2)
        B(2,10) = dNdx(3,2)
        B(2,14) = dNdx(4,2)
        B(2,18) = dNdx(5,2)
        B(2,20) = dNdx(6,2)
        B(2,22) = dNdx(7,2)
        B(2,24) = dNdx(8,2)

        B(3,1)  = dNdx(1,2)
        B(3,5)  = dNdx(2,2)
        B(3,9)  = dNdx(3,2)
        B(3,13) = dNdx(4,2)
        B(3,17) = dNdx(5,2)
        B(3,19) = dNdx(6,2)
        B(3,21) = dNdx(7,2)
        B(3,23) = dNdx(8,2)

        B(4,2)  = dNdx(1,1)
        B(4,6)  = dNdx(2,1)
        B(4,10) = dNdx(3,1)
        B(4,14) = dNdx(4,1)
        B(4,18) = dNdx(5,1)
        B(4,20) = dNdx(6,1)
        B(4,22) = dNdx(7,1)
        B(4,24) = dNdx(8,1)
        
        B(4,3)  = N_ph(1)
        B(4,7)  = N_ph(2)
        B(4,11) = N_ph(3)
        B(4,15) = N_ph(4)

        B(5,4)  = N_ph(1)
        B(5,8)  = N_ph(2)
        B(5,12) = N_ph(3)
        B(5,16) = N_ph(4)
        
        B(6,3)  = dNdx_ph(1, 1)
        B(6,7)  = dNdx_ph(2, 1)
        B(6,11) = dNdx_ph(3, 1)
        B(6,15) = dNdx_ph(4, 1)

        B(7,3)  = dNdx_ph(1, 2)
        B(7,7)  = dNdx_ph(2, 2)
        B(7,11) = dNdx_ph(3, 2)
        B(7,15) = dNdx_ph(4, 2)

        B(8,4)  = dNdx_ph(1, 1)
        B(8,8)  = dNdx_ph(2, 1)
        B(8,12) = dNdx_ph(3, 1)
        B(8,16) = dNdx_ph(4, 1)

        B(9,4)  = dNdx_ph(1, 2)
        B(9,8)  = dNdx_ph(2, 2)
        B(9,12) = dNdx_ph(3, 2)
        B(9,16) = dNdx_ph(4, 2)

        strain = matmul(B,U)
        dstrain = matmul(B,DU(1:2*NNODE+2*4,1))

        con = strain(5)

        elstrain(1) = strain(1)
        elstrain(2) = strain(2)
        elstrain(3) = 0.D0
        elstrain(4) = strain(3)
        elstress = matmul(Del,(elstrain - con*Om*delta/3.D0))

        q(1) = elstress(1)
        q(2) = elstress(2)
        q(3) = elstress(4)
        q(4) = strain(4) - 2.D0*WC*con*(con - 1.D0)*(2.D0*con - 1.D0)
     1    - Om*(elstress(1)+elstress(2)+elstress(3))/3.D0

        q(5) = dstrain(5)/DTIME
        q(6) = -kappa*strain(8)
        q(7) = -kappa*strain(9)
        q(8) = dc*(strain(6) + (theta - 1.D0)*dstrain(6))
        q(9) = dc*(strain(7) + (theta - 1.D0)*dstrain(7))

!       print*, D

        D(4,5) = -WC*(12.D0*con*con - 12.D0*con + 2.D0)
     1    + Om*Om*(Del(1,1)+Del(1,2)+Del(1,3)+Del(2,1)+Del(2,2)+Del(2,3)
     2    + Del(3,1)+Del(3,2)+Del(3,3))/9.D0     

        RHS(1:2*NNODE+2*4,1) = RHS(1:2*NNODE+2*4,1)
     1     - matmul(transpose(B(1:9,1:2*NNODE+2*4)),q)*
     2                                            w(kint)*determinant

     
        AMATRX(1:2*NNODE+2*4,1:2*NNODE+2*4) = 
     1   AMATRX(1:2*NNODE+2*4,1:2*NNODE+2*4) +
     2   matmul(transpose(B(1:9,1:2*NNODE+2*4))
     3   , matmul(D,B(1:9,1:2*NNODE+2*4)))*w(kint)*determinant

        ENERGY(2) = ENERGY(2)
     1     + 0.5D0*dot_product(elstress,elstrain)*w(kint)*determinant           ! Store the elastic strain energy
        
       if (NSVARS>=n_points*4) then   ! Store stress at each integration point (if space was allocated to do so)
          SVARS(4*kint-3:4*kint) = elstress(1:4)
        endif
      end do
      PNEWDT = 1.d0          ! This leaves the timestep unchanged (ABAQUS will use its own algorithm to determine DTIME)
    
      elseif (JTYPE==2) then

      print*, 'Incompatible mode not implemented'

    !   print*, 'The incompatible mode is running'


    !   call abq_UEL_2D_integrationpoints(1, NNODE, xm, wm)

    !   call abq_UEL_2D_shapefunctions(xm,NNODE,N,dNdxi)
    !   dxdxi = matmul(coords(1:2,1:NNODE),dNdxi(1:NNODE,1:2))
    !   call abq_UEL_invert2d(dxdxi,dxidx,det0)

    !   call abq_UEL_2D_integrationpoints(n_points, NNODE, xi, w)

    !   do kint = 1, n_points
    !     call abq_UEL_2D_shapefunctions(xi(1:2,kint),NNODE,N,dNdxi)
    !     dxdxi = matmul(coords(1:2,1:NNODE),dNdxi(1:NNODE,1:2))
    !     call abq_UEL_invert2d(dxdxi,dxidx,determinant)
    !     dNdx(1:NNODE,1:2) = matmul(dNdxi(1:NNODE,1:2),dxidx)

    !     B = 0.d0
    !     B(1,1:2*NNODE-1:2) = dNdx(1:NNODE,1)
    !     B(2,2:2*NNODE:2) = dNdx(1:NNODE,2)
    !     B(4,1:2*NNODE-1:2) = dNdx(1:NNODE,2)
    !     B(4,2:2*NNODE:2) = dNdx(1:NNODE,1)

    !     B(1,2*NNODE+1) = xi(1,kint)*(det0/determinant)*dxidx(1,1)
    !     B(1,2*NNODE+3) = xi(2,kint)*(det0/determinant)*dxidx(2,1)
    !     B(2,2*NNODE+2) = xi(1,kint)*(det0/determinant)*dxidx(1,2)
    !     B(2,2*NNODE+4) = xi(2,kint)*(det0/determinant)*dxidx(2,2)
    !     B(4,2*NNODE+1) = B(2,(2*NNODE)+2)
    !     B(4,2*NNODE+2) = B(1,(2*NNODE)+1)
    !     B(4,2*NNODE+3) = B(2,(2*NNODE)+4)
    !     B(4,2*NNODE+4) = B(1,(2*NNODE)+3)

    !     ktemp(1:2*NNODE+4,1:2*NNODE+4)=ktemp(1:2*NNODE+4,1:(2*NNODE)+4) 
    !  1   + matmul(transpose(B(1:4,1:2*NNODE+4)),
    !  2     matmul(D,B(1:4,1:2*NNODE+4)))*w(kint)*determinant

    !   enddo 
              
    !   kaa(1:4,1:4) = ktemp(2*NNODE+1:2*NNODE+4,2*NNODE+1:2*NNODE+4)
    !   kuu(1:2*NNODE,1:2*NNODE) = ktemp(1:2*NNODE,1:2*NNODE)        
    !   kau(1:4,1:2*NNODE) = ktemp(2*NNODE+1:2*NNODE+4,1:2*NNODE)
    !   kua(1:2*NNODE,1:4) = ktemp(1:2*NNODE,2*NNODE+1:2*NNODE+4)

    !   call abq_inverse_LU(kaa,kaainv,4)
      
    !   alpha = -(matmul(kaainv,matmul(kau(1:4,1:2*NNODE),U(1:2*NNODE))))
      
    !   do kint = 1, n_points
    !     call abq_UEL_2D_shapefunctions(xi(1:2,kint),NNODE,N,dNdxi)
    !     dxdxi = matmul(coords(1:2,1:NNODE),dNdxi(1:NNODE,1:2))
    !     call abq_UEL_invert2d(dxdxi,dxidx,determinant)
    !     dNdx(1:NNODE,1:2) = matmul(dNdxi(1:NNODE,1:2),dxidx)

    !     B = 0.d0
    !     B(1,1:2*NNODE-1:2) = dNdx(1:NNODE,1)
    !     B(2,2:2*NNODE:2) = dNdx(1:NNODE,2)
    !     B(4,1:2*NNODE-1:2) = dNdx(1:NNODE,2)
    !     B(4,2:2*NNODE:2) = dNdx(1:NNODE,1)

    !     B(1,(2*NNODE)+1) = xi(1,kint)*(det0/determinant)*dxidx(1,1)
    !     B(1,(2*NNODE)+3) = xi(2,kint)*(det0/determinant)*dxidx(2,1)
    !     B(2,(2*NNODE)+2) = xi(1,kint)*(det0/determinant)*dxidx(1,2)
    !     B(2,(2*NNODE)+4) = xi(2,kint)*(det0/determinant)*dxidx(2,2)
    !     B(4,(2*NNODE)+1) = B(2,(2*NNODE)+2)
    !     B(4,(2*NNODE)+2) = B(1,(2*NNODE)+1)
    !     B(4,(2*NNODE)+3) = B(2,(2*NNODE)+4)
    !     B(4,(2*NNODE)+4) = B(1,(2*NNODE)+3)

    !     strain =matmul(B(1:4,1:2*NNODE+4),[U(1:2*NNODE),alpha(1:4)])
    !     stress = matmul(D,strain)

    !     rhs_temp(1:2*NNODE+4) = rhs_temp(1:2*NNODE+4)
    !  1    - matmul(transpose(B(1:4,1:2*NNODE+4)),stress(1:4))*
    !  2                                          w(kint)*determinant
      
    !     ENERGY(2) = ENERGY(2)
    !  1     + 0.5D0*dot_product(stress,strain)*w(kint)*determinant           ! Store the elastic strain energy

    !     if (NSVARS>=n_points*4) then   ! Store stress at each integration point (if space was allocated to do so)
    !       SVARS(4*kint-3:4*kint) = stress(1:4)
    !     endif

    !   enddo 

    !   AMATRX(1:2*NNODE,1:2*NNODE) = kuu(1:2*NNODE,1:2*NNODE)
    !  1     -matmul(kua(1:2*NNODE,1:4),matmul(kaainv,kau(1:4,1:2*NNODE)))

    !   RHS(1:2*NNODE,1)= rhs_temp(1:2*NNODE)-matmul(kua(1:2*NNODE,1:4),
    !  1             matmul(kaainv,rhs_temp(2*NNODE+1:2*NNODE+4)))
        
       endif
    
    !
    !   Apply distributed loads
    !
    !   Distributed loads are specified in the input file using the Un option in the input file.
    !   n specifies the face number, following the ABAQUS convention
    !
  !**************************

    !   do j = 1,NDLOAD

    !     call abq_facenodes_2D(NNODE,iabs(JDLTYP(j,1)),
    !  1                                     face_node_list,nfacenodes)
    !     do i = 1,nfacenodes
    !         face_coords(1:2,i) = coords(1:2,face_node_list(i))
    !     end do

    !     if (nfacenodes == 2) n_points = 2
    !     if (nfacenodes == 3) n_points = 3

    !     call abq_UEL_1D_integrationpoints(n_points, nfacenodes, xi1, w1)
        
    !     do kint = 1,n_points
    !       call abq_UEL_1D_shapefunctions(xi1(kint),
    !  1                        nfacenodes,N1,dNdxi1)
    !       dxdxi1 = matmul(face_coords(1:2,1:nfacenodes),
    !  1            dNdxi1(1:nfacenodes))

    !       norm(1)=dxdxi1(2)
    !       norm(2)=-dxdxi1(1)
          
    !       do i = 1,nfacenodes
    !         ipoin = 2*face_node_list(i)-1
    !         RHS(ipoin:ipoin+1,1) = RHS(ipoin:ipoin+1,1)
    !  1        - N1(1:nfacenodes)*adlmag(j,1)*norm(1:2)*w1(kint)  ! Note determinant is already in normal
                
    !       end do
    !     end do
    !   end do
        
      return

      END SUBROUTINE UEL_PH

!      INCLUDE 'facenodes.for'
!      INCLUDE 'shapefunctions.for'
!      INCLUDE 'integrationpoints.for'
!      INCLUDE 'elementutilities.for'