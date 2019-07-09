!
!    ABAQUS format UEL subroutine
!
!    This file is compatible with both EN234_FEA and ABAQUS/Standard
!
!    The example implements a 2D Continuum Beam element
!=========================== ABAQUS format user element subroutine ===================

      SUBROUTINE UEL_CB(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
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
      integer      :: i,j,n_points,kint
    !
      double precision  ::  xi(2,5)                          ! Area integration points
      double precision  ::  w(9)                             ! Area integration weights
      double precision  ::  N(9)                             ! 2D shape functions
      double precision  ::  dNdxi(9,2)                       ! 2D shape function derivatives
      double precision  ::  dNdx(9,2)                        ! Spatial derivatives
      double precision  ::  dxdxi(2,2)                       ! Derivative of spatial coords wrt normalized coords
    !
      double precision  ::  strain(3)                         ! Strain vector contains [e11, e22, 2e12]
      double precision  ::  strain_hat(2)                     ! strain hat vector contains [ehat11, 2ehat12]
      double precision  ::  ehat(2,2)                         ! ehat for constructing R matrix
      double precision  ::  stress(2)                         ! Stress vector contains [s11, s22]
      double precision  ::  D(9,9)                            ! stress = D*(strain)  (NOTE FACTOR OF 2 in shear strain)
      double precision  ::  B(9,22)                           ! strain = B*(dof_total)
      double precision  ::  dxidx(2,2), determinant           ! Jacobian inverse and determinant
      double precision  ::  E, xnu, xh, xw                    ! Material properties
      double precision  ::  L                                 ! 

      double precision  ::  sine, cosine                      ! 
      double precision  ::  x(2,6)                            ! 
      double precision  ::  T(8,6), R(2,3)                    ! Rotation and Transformation tensor
      double precision  ::  axial, shear, moment


    !
    !     Example ABAQUS UEL implementing 2D linear elastic elements
    !     Includes option for incompatible mode elements
    !     El props are:

    !     PROPS(1)         Beams depth
    !     PROPS(2)         Beams width
    !     PROPS(3)         Young's modulus
    !     PROPS(4)         Poisson's ratio

    

      xh = PROPS(1)
      xw = PROPS(2)
      E = PROPS(3)
      xnu = PROPS(4)

      RHS(1:MLVARX,1) = 0.d0
      AMATRX(1:NDOFEL,1:NDOFEL) = 0.d0

      SVARS(1:NSVARS) = 0.d0

      axial = 0.d0
      shear = 0.d0
      moment = 0.d0

      xi = 0.D0
      n_points = 5

      call abq_UEL_1D_integrationpoints(n_points, 4, xi(2,1:5), w)

      x = 0.D0
    
      ! Master Coordinates

      x(1:2, 5) = COORDS(1:2,1)
      x(1:2, 6) = COORDS(1:2,2)

      L = sqrt((x(1,5)-x(1,6))**2.D0+(x(2,5)-x(2,6))**2.D0)

      sine = (x(2,6)-x(2,5))/L
      cosine = (x(1,6)-x(1,5))/L

      ! Slave Coordinates

      x(1,1) = x(1,5)+(xh/2.d0)*sine
      x(2,1) = x(2,5)-(xh/2.d0)*cosine
      x(1,2) = x(1,6)+(xh/2.d0)*sine
      x(2,2) = x(2,6)-(xh/2.d0)*cosine
      x(1,3) = x(1,6)-(xh/2.d0)*sine
      x(2,3) = x(2,6)+(xh/2.d0)*cosine
      x(1,4) = x(1,5)-(xh/2.d0)*sine
      x(2,4) = x(2,5)+(xh/2.d0)*cosine

      ! T Computation

      T=0.d0
      T(1,1)=1.d0
      T(1,3)=x(2,5)-x(2,1)
      T(2,2)=1.d0
      T(2,3)= -1.d0*(x(1,5)-x(1,1))
      T(3,4)=1.d0
      T(3,6)=x(2,6)-x(2,2)
      T(4,5)=1.d0
      T(4,6)=x(1,6)-x(1,2)
      T(5,4)=1.d0
      T(5,6)=x(2,6)-x(2,3)
      T(6,5)=1.d0
      T(6,6)=x(1,6)-x(1,3)
      T(7,1)=1.d0
      T(7,3)=x(2,5)-x(2,4)
      T(8,2)=1.d0
      T(8,3)= -1.d0*(x(1,5)-x(1,4))
    
      ! D computation

      D = 0.d0
      D(1,1)=E
      D(2,2)=0.5d0*E/(1+xnu)

      ENERGY(1:8) = 0.d0
    
!      JTYPE = 2

      if (JTYPE==3) then

      print*, 'The Continuum Beam Element'

      do kint = 1, n_points
        call abq_UEL_2D_shapefunctions(xi(1:2,kint),4,N,dNdxi)
        dxdxi = matmul(x(1:2,1:4),dNdxi(1:4,1:2))
        call abq_UEL_invert2d(dxdxi,dxidx,determinant)
        dNdx(1:4,1:2) = matmul(dNdxi(1:4,1:2),dxidx)


        B = 0.d0
        B(1,1:2*NNODE:2)=dNdx(1:4,1)
        B(2,2:2*NNODE:2)=dNdx(1:4,2)
        B(3,1:2*NNODE:2)=dNdx(1:4,2)
        B(3,2:2*NNODE:2)=dNdx(1:4,1)

        !R Computation
        ehat(1:2,1)=dxdxi(1:2,1)/norm2(dxdxi(1:2,1))
        ehat(1:2,2)=dxdxi(1:2,2)/norm2(dxdxi(1:2,2))

        R = 0.d0
        R(1,1)=ehat(1,1)**2
        R(1,2)=ehat(1,2)**2
        R(1,3)=ehat(1,1)*ehat(1,2)
        R(2,1)=2.d0*ehat(1,1)*ehat(2,1)
        R(2,2)=2.d0*ehat(1,2)*ehat(2,2)
        R(2,3)=ehat(1,1)*ehat(2,2)+ehat(1,2)*ehat(2,1)

        strain=matmul(B(1:3,1:8),matmul(T(1:8,1:6),U(1:6)))
        strain_hat(1:2)=matmul(R(1:2,1:3),strain(1:3))
        stress(1:2)=matmul(D(1:2,1:2),strain_hat(1:2))


        RHS(1:6,1) = RHS(1:6,1)
     1   -(L*xh*xw*0.5d0)*matmul(transpose(matmul(R(1:2,1:3),
     2   matmul(B(1:3,1:8),T(1:8,1:6)))),stress(1:2))*w(kint)

        AMATRX(1:6,1:6) = AMATRX(1:6,1:6)
     1   + (L*xh*xw*0.5d0)*matmul(transpose(matmul(R(1:2,1:3),
     2   matmul(B(1:3,1:8),T(1:8,1:6)))),matmul(D(1:2,1:2),
     3   matmul(R(1:2,1:3),matmul(B(1:3,1:8),T(1:8,1:6)))))*w(kint)
	   
        axial = axial+(xw*xh*0.5d0)*stress(1)*w(kint)
        shear = shear+(xw*xh*0.5d0)*stress(2)*w(kint)
        moment = moment+(xw*xh**2*0.25d0)*stress(1)*xi(2,kint)*w(kint)

!        ENERGY(2) = ENERGY(2)
!     1     + 0.5D0*dot_product(stress,strain)*w(kint)*determinant           ! Store the elastic strain energy

        SVARS(1) = axial
        SVARS(2) = shear
        SVARS(3) = moment
      end do

      PNEWDT = 1.d0          ! This leaves the timestep unchanged (ABAQUS will use its own algorithm to determine DTIME)
      else
        print*, 'Wrong element is specified.'
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

      END SUBROUTINE UEL_CB


!      INCLUDE 'facenodes.for'
!      INCLUDE 'shapefunctions.for'
!      INCLUDE 'integrationpoints.for'
!      INCLUDE 'elementutilities.for'