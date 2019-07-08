!
!    ABAQUS format UEL subroutine
!
!    This file is compatible with both EN234_FEA and ABAQUS/Standard
!
!    The example implements a standard fully integrated 3D linear elastic continuum element
!
!    The file also contains the following subrouines:
!          abq_UEL_invert3D                       - computes the inverse and determinant of a 3x3 matrix
!
!=========================== ABAQUS format user element subroutine ===================

      SUBROUTINE UEL_3D(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
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
      integer      :: i,j,n_points,kint, nfacenodes, ipoin
      integer      :: face_node_list(8)                       ! List of nodes on an element face
    !
      double precision  ::  xi(3,64)                          ! Volumetric Integration points
      double precision  ::  w(64)                             ! Integration weights
      double precision  ::  N(20)                             ! 3D Shape functions
      double precision  ::  dNdxi(20,3)                       ! 3D Shape function derivatives
      double precision  ::  dxdxi(3,3)                        ! Derivative of position wrt normalized coords
      double precision  ::  dNdx(20,3)                        ! Derivative of shape functions wrt spatial coords
    !
    !   Variables below are for computing integrals over element faces
      double precision  ::  face_coords(3,8)                  ! Coords of nodes on an element face
      double precision  ::  xi2(2,9)                          ! Area integration points
      double precision  ::  N2(9)                             ! 2D shape functions
      double precision  ::  dNdxi2(9,2)                       ! 2D shape function derivatives
      double precision  ::  norm(3)                           ! Normal to an element face
      double precision  ::  dxdxi2(3,2)                       ! Derivative of spatial coord wrt normalized areal coord
    !
      double precision  ::  strain(6)                         ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
      double precision  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
      double precision  ::  D(6,6)                            ! stress = D*(strain)  (NOTE FACTOR OF 2 in shear strain)
      double precision  ::  B(6,60)                           ! strain = B*(dof_total)
      double precision  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
      double precision  ::  E, xnu, D44, D11, D12             ! Material properties

    !
    !     Example ABAQUS UEL implementing 3D linear elastic elements
    !     El props are:

    !     PROPS(1)         Young's modulus
    !     PROPS(2)         Poisson's ratio

      if (NNODE == 4) n_points = 1               ! Linear tet
      if (NNODE == 10) n_points = 4              ! Quadratic tet
      if (NNODE == 8) n_points = 8               ! Linear Hex
      if (NNODE == 20) n_points = 27             ! Quadratic hex

      call abq_UEL_3D_integrationpoints(n_points, NNODE, xi, w)

      if (MLVARX<3*NNODE) then
        write(6,*) ' Error in abaqus UEL '
        write(6,*) ' Variable MLVARX must exceed 3*NNODE'
        write(6,*) ' MLVARX = ',MLVARX,' NNODE = ',NNODE
        stop
      endif

      RHS(1:MLVARX,1) = 0.d0
      AMATRX(1:NDOFEL,1:NDOFEL) = 0.d0

      D = 0.d0
      E = PROPS(1)
      xnu = PROPS(2)
      d44 = 0.5D0*E/(1+xnu)
      d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
      d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
      D(1:3,1:3) = d12
      D(1,1) = d11
      D(2,2) = d11
      D(3,3) = d11
      D(4,4) = d44
      D(5,5) = d44
      D(6,6) = d44

      ENERGY(1:8) = 0.d0

    !     --  Loop over integration points
      do kint = 1, n_points
        call abq_UEL_3D_shapefunctions(xi(1:3,kint),NNODE,N,dNdxi)
        dxdxi = matmul(coords(1:3,1:NNODE),dNdxi(1:NNODE,1:3))
        call abq_UEL_invert3d(dxdxi,dxidx,determinant)
        dNdx(1:NNODE,1:3) = matmul(dNdxi(1:NNODE,1:3),dxidx)
        B = 0.d0
        B(1,1:3*NNODE-2:3) = dNdx(1:NNODE,1)
        B(2,2:3*NNODE-1:3) = dNdx(1:NNODE,2)
        B(3,3:3*NNODE:3)   = dNdx(1:NNODE,3)
        B(4,1:3*NNODE-2:3) = dNdx(1:NNODE,2)
        B(4,2:3*NNODE-1:3) = dNdx(1:NNODE,1)
        B(5,1:3*NNODE-2:3) = dNdx(1:NNODE,3)
        B(5,3:3*NNODE:3)   = dNdx(1:NNODE,1)
        B(6,2:3*NNODE-1:3) = dNdx(1:NNODE,3)
        B(6,3:3*NNODE:3)   = dNdx(1:NNODE,2)

        strain = matmul(B(1:6,1:3*NNODE),U(1:3*NNODE))

        stress = matmul(D,strain)
        RHS(1:3*NNODE,1) = RHS(1:3*NNODE,1)
     1   - matmul(transpose(B(1:6,1:3*NNODE)),stress(1:6))*
     2                                          w(kint)*determinant

        AMATRX(1:3*NNODE,1:3*NNODE) = AMATRX(1:3*NNODE,1:3*NNODE)
     1  + matmul(transpose(B(1:6,1:3*NNODE)),matmul(D,B(1:6,1:3*NNODE)))
     2                                             *w(kint)*determinant

        ENERGY(2) = ENERGY(2)
     1   + 0.5D0*dot_product(stress,strain)*w(kint)*determinant           ! Store the elastic strain energy

        if (NSVARS>=n_points*6) then   ! Store stress at each integration point (if space was allocated to do so)
            SVARS(6*kint-5:6*kint) = stress(1:6)
        endif
      end do


      PNEWDT = 1.d0          ! This leaves the timestep unchanged (ABAQUS will use its own algorithm to determine DTIME)
    !
    !   Apply distributed loads
    !
    !   Distributed loads are specified in the input file using the Un option in the input file.
    !   n specifies the face number, following the ABAQUS convention
    !
    !
      do j = 1,NDLOAD

        call abq_facenodes_3D(NNODE,iabs(JDLTYP(j,1)),
     1                                     face_node_list,nfacenodes)

        do i = 1,nfacenodes
            face_coords(1:3,i) = coords(1:3,face_node_list(i))
        end do

        if (nfacenodes == 3) n_points = 3
        if (nfacenodes == 6) n_points = 4
        if (nfacenodes == 4) n_points = 4
        if (nfacenodes == 8) n_points = 9

        call abq_UEL_2D_integrationpoints(n_points, nfacenodes, xi2, w)

        do kint = 1,n_points
            call abq_UEL_2D_shapefunctions(xi2(1:2,kint),
     1                        nfacenodes,N2,dNdxi2)
            dxdxi2 = matmul(face_coords(1:3,1:nfacenodes),
     1                           dNdxi2(1:nfacenodes,1:2))
            norm(1)=(dxdxi2(2,1)*dxdxi2(3,2))-(dxdxi2(2,2)*dxdxi2(3,1))
            norm(2)=(dxdxi2(1,1)*dxdxi2(3,2))-(dxdxi2(1,2)*dxdxi2(3,1))
            norm(3)=(dxdxi2(1,1)*dxdxi2(2,2))-(dxdxi2(1,2)*dxdxi2(2,1))

            do i = 1,nfacenodes
                ipoin = 3*face_node_list(i)-2
                RHS(ipoin:ipoin+2,1) = RHS(ipoin:ipoin+2,1)
     1                 - N2(1:nfacenodes)*adlmag(j,1)*norm(1:3)*w(kint)      ! Note determinant is already in normal
            end do
        end do
      end do
      return

      END SUBROUTINE UEL_3D

!      INCLUDE 'facenodes.for'
!      INCLUDE 'shapefunctions.for'
!      INCLUDE 'integrationpoints.for'
!      INCLUDE 'elementutilities.for'
