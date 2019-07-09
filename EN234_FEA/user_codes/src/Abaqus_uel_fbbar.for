!
!    ABAQUS format UEL subroutine
!
!    This file is compatible with both EN234_FEA and ABAQUS/Standard
!
!    3D finite-strain B-Bar element with a Neo-Hookean material model
!
!    The file also contains the following subroutines: 
!          abq_UEL_3D_integrationpoints           - defines integration points for 3D continuum elements
!          abq_UEL_3D_shapefunctions              - defines shape functions for 3D continuum elements
!          abq_UEL_invert3D                       - computes the inverse and determinant of a 3x3 matrix
!          abq_facenodes_3D                       - returns list of nodes on the face of a 3D element
!
!=========================== ABAQUS format user element subroutine ===================

      SUBROUTINE UEL_BB(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
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
      double precision  ::  dNdy(20,3)                        ! Derivative of shape functions wrt deformed coords
      double precision  ::  dNdyvec(60)                       ! Shape fct derivs in vector form
      double precision  ::  dNbardy(20,3)                     ! Volume-averaged shape function derivatives
      double precision  ::  dNbardyvec(60)                    ! Volume-averaged shape function derivatives
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
      double precision  ::  skk                               ! Trace of Kirchoff stress
      double precision  ::  delta(3,3), deltavec(6)           ! Identity matrix, and its vector form
      double precision  ::  uint(3,20)                        ! Internal displacement matrix for calculating F via matmul (DOF,node)
      double precision  ::  F(3,3), Finv(3,3)                 ! Deformation gradient and its inverse
      double precision  ::  Fbar(3,3)                         ! Vol avgd deformation gradient
      double precision  ::  Ja                                ! Determinant of deformation gradient
      double precision  ::  Jbar                              ! Volume-averaged determinant
      double precision  ::  Vel                               ! Element volume
      double precision  ::  Sigma(60,60)                      ! Geometric stiffness
      double precision  ::  Q(60,60)                          ! Geometric stiffness
      double precision  ::  P(60,60)                          ! Geometric stiffness
      double precision  ::  Pvec(60)                          ! Vector used in calculating P
      double precision  ::  Pmat(60,60)                       ! Matrix used in calculating P  
      double precision  ::  LCG(3,3), LCGvec(6), LCGkk        ! Left cauchy green stretch tensor, its vector form, and its trace
      double precision  ::  LCGi(3,3), LCGivec(6), LCGdet     ! Inverses of LCG tensor
      double precision  ::  Kirch(3,3)                        ! Kirchoff stress tensor
      double precision  ::  D(6,6)                            ! Tangent stiffness matrix
      double precision  ::  Bbar(6,60), Bbar2(6,60)           ! Modified B matrix       
      double precision  ::  Bstar(9,60), Bstar2(9,60)         ! Second modified B matrix
      double precision  ::  G(6,9)                            ! Used in calculating stiffness
      double precision  ::  S(3,20), Svec(60), Smat(60,60)    ! Used for computing sigma
      double precision  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
      double precision  ::  mu, K                             ! Material properties
      
    !
    ! ABAQUS UEL implementing 3D finite-strain B-bar element with a NeoHookian material model
    !     
      
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
      
      ! Initialize variables
      mu = PROPS(1)
      K  = PROPS(2)
      
      delta       = 0.d0
      deltavec    = 0.d0
      delta(1,1)  = 1.d0
      delta(2,2)  = 1.d0
      delta(3,3)  = 1.d0
      deltavec(1) = 1.d0
      deltavec(2) = 1.d0
      deltavec(3) = 1.d0
      
      uint = 0.d0
      do j = 0, (NNODE-1)
          do i = 1, 3
              uint(i,j+1) = U(3*j + i)
          end do
      end do

      ENERGY(1:8) = 0.d0

    !!!First loop over integration points - determine volume-averaged functions
      dNbardy = 0.d0
      Jbar    = 0.d0
      Vel     = 0.d0
      P       = 0.d0
      
      do kint = 1, n_points
        
        call abq_UEL_3D_shapefunctions(xi(1:3,kint),NNODE,N,dNdxi)
        dxdxi = matmul(COORDS(1:3,1:NNODE),dNdxi(1:NNODE,1:3))
        
        call abq_UEL_invert3d(dxdxi,dxidx,determinant)
        dNdx(1:NNODE,1:3) = matmul(dNdxi(1:NNODE,1:3),dxidx)
        
        F    = 0.d0
        Finv = 0.d0
        Ja   = 0.d0
        F = matmul(uint(1:3,1:NNODE),dNdx(1:NNODE,1:3))+delta
        call abq_UEL_invert3d(F,Finv,Ja)
        
        dNdy = 0.d0
        dNdy(1:NNODE,1:3) = matmul(dNdx(1:NNODE,1:3),Finv)
        
        !Add contributions from each integration point to following terms:
        dNbardy(1:NNODE,1:3) = dNbardy(1:NNODE,1:3) 
     1                        + Ja*dNdy(1:NNODE,1:3)*w(kint)*determinant
        Jbar                 = Jbar + Ja*w(kint)*determinant
        Vel                  = Vel  + w(kint)*determinant
        
        !Add contributions from each integration point to integral in P
        dNdyvec = 0.d0
        Pvec = 0.d0
        Pmat = 0.d0
        dNdyvec(1:3*NNODE) = reshape(transpose(dNdy(1:NNODE,1:3)),
     1                                                    (/3*NNODE/))
        do i = 1,NNODE
            Pvec(1:3*NNODE) = reshape(spread(transpose(dNdy(i:i,1:3)),
     1       dim=2,ncopies=NNODE),(/3*NNODE/))
            Pmat(3*i-2:3*i,1:3*NNODE) = spread(Pvec,dim=1,ncopies=3)
        end do
        
        P = P + (Ja*spread(dNdyvec,dim=2,ncopies=3*NNODE)
     1  *spread(dNdyvec,dim=1,ncopies=3*NNODE) -Ja*Pmat*transpose(Pmat))
     2                                              *w(kint)*determinant
        
      end do
      
      !Normalize by Jbar*Vel
      Jbar    = Jbar/Vel   
      dNbardy = dNbardy/(Jbar*Vel)
      P       = P/(Jbar*Vel)
      
      !Add last term to P
      dNbardyvec = 0.d0
      dNbardyvec(1:3*NNODE) = reshape(transpose(dNbardy(1:NNODE,1:3)),
     1 (/3*NNODE/))
      
      P = P - (spread(dNbardyvec,dim=2,ncopies=3*NNODE)*
     1         spread(dNbardyvec,dim=1,ncopies=3*NNODE))
      
    !!!Second loop over integration points - stiffness, rhs
      N     = 0.d0
      dNdxi = 0.d0
      dxdxi = 0.d0
      dxidx = 0.d0
      dNdx  = 0.d0
      
      do kint = 1, n_points
        
        call abq_UEL_3D_shapefunctions(xi(1:3,kint),NNODE,N,dNdxi)
        dxdxi = matmul(COORDS(1:3,1:NNODE),dNdxi(1:NNODE,1:3))
        
        call abq_UEL_invert3d(dxdxi,dxidx,determinant)
        dNdx(1:NNODE,1:3) = matmul(dNdxi(1:NNODE,1:3),dxidx)
        
        F    = 0.d0
        Finv = 0.d0
        Ja   = 0.d0
        F = matmul(uint(1:3,1:NNODE),dNdx(1:NNODE,1:3))+delta
        call abq_UEL_invert3d(F,Finv,Ja)
        
        dNdy = 0.d0
        dNdy(1:NNODE,1:3) = matmul(dNdx(1:NNODE,1:3),Finv)
        
        Fbar = 0.d0
        Fbar = ((Jbar/Ja)**(1.d0/3.d0))*F
        
        !Assemble kirchoff stress, tangent stiffness matrix
        LCG     = 0.d0
        LCGi    = 0.d0
        LCGvec  = 0.d0
        LCGivec = 0.d0
        LCGkk   = 0.d0
        LCG     = matmul(Fbar,transpose(Fbar))
        call abq_UEL_invert3d(LCG,LCGi,LCGdet)
        
        LCGvec(1) = LCG(1,1)
        LCGvec(2) = LCG(2,2)
        LCGvec(3) = LCG(3,3)
        LCGvec(4) = LCG(1,2)
        LCGvec(5) = LCG(1,3)
        LCGvec(6) = LCG(2,3)
        LCGivec(1) = LCGi(1,1)
        LCGivec(2) = LCGi(2,2)
        LCGivec(3) = LCGi(3,3)
        LCGivec(4) = LCGi(1,2)
        LCGivec(5) = LCGi(1,3)
        LCGivec(6) = LCGi(2,3)
        LCGkk = LCG(1,1) + LCG(2,2) + LCG(3,3)
        
        Kirch = 0.d0
        Kirch = (mu/(Jbar**(2.d0/3.d0)))*(LCG - (LCGkk/3.d0)*delta)
     1           + K*Jbar*(Jbar - 1.d0)*delta
        stress = 0.d0
        stress(1) = Kirch(1,1)
        stress(2) = Kirch(2,2)
        stress(3) = Kirch(3,3)
        stress(4) = Kirch(1,2)
        stress(5) = Kirch(1,3)
        stress(6) = Kirch(2,3)
        
        D = 0.d0
        D(1:3,1:3) = delta
        D(4:6,4:6) = 0.5d0*delta
        D = (mu/(Jbar**(2.d0/3.d0)))*D
        D = D + 
     &  (mu/(3.d0*Jbar**(2.d0/3.d0)))*((LCGkk/3.d0)*
     &  spread(deltavec,dim=2,ncopies=6)*spread(LCGivec,dim=1,ncopies=6)
     &  -spread(deltavec,dim=2,ncopies=6)*
     &   spread(deltavec,dim=1,ncopies=6)
     &  -spread(LCGvec,dim=2,ncopies=6)*spread(LCGivec,dim=1,ncopies=6))
     &  + K*Jbar*(Jbar-0.5d0)*spread(deltavec,dim=2,ncopies=6)*
     &   spread(LCGivec,dim=1,ncopies=6)
        
        skk = 0.d0
        skk = stress(1) + stress(2) + stress(3)
        
        
        !Assemble G
        G = 0.d0
        G(1,1) = 2.d0*LCG(1,1)
        G(1,4) = 2.d0*LCG(1,2)
        G(1,6) = 2.d0*LCG(1,3)
        G(2,2) = 2.d0*LCG(2,2)
        G(2,5) = 2.d0*LCG(1,2)
        G(2,8) = 2.d0*LCG(2,3)
        G(3,3) = 2.d0*LCG(3,3)
        G(3,7) = 2.d0*LCG(1,3)
        G(3,9) = 2.d0*LCG(2,3)
        G(4,1) = 2.d0*LCG(1,2)
        G(4,2) = 2.d0*LCG(1,2)
        G(4,4) = 2.d0*LCG(2,2)
        G(4,5) = 2.d0*LCG(1,1)
        G(4,6) = 2.d0*LCG(2,3)
        G(4,8) = 2.d0*LCG(1,3)
        G(5,1) = 2.d0*LCG(1,3)
        G(5,3) = 2.d0*LCG(1,3)
        G(5,4) = 2.d0*LCG(2,3)
        G(5,6) = 2.d0*LCG(3,3)
        G(5,7) = 2.d0*LCG(1,1)
        G(5,9) = 2.d0*LCG(1,2)
        G(6,2) = 2.d0*LCG(2,3)
        G(6,3) = 2.d0*LCG(2,3)
        G(6,5) = 2.d0*LCG(1,3)
        G(6,7) = 2.d0*LCG(1,2)
        G(6,8) = 2.d0*LCG(3,3)
        G(6,9) = 2.d0*LCG(2,2)
        
        !Assemble Bbar
        Bbar = 0.d0
        Bbar(1,1:3*NNODE-2:3) = dNdy(1:NNODE,1)
        Bbar(2,2:3*NNODE-1:3) = dNdy(1:NNODE,2)
        Bbar(3,3:3*NNODE:3)   = dNdy(1:NNODE,3)
        Bbar(4,1:3*NNODE-2:3) = dNdy(1:NNODE,2)
        Bbar(4,2:3*NNODE-1:3) = dNdy(1:NNODE,1)
        Bbar(5,1:3*NNODE-2:3) = dNdy(1:NNODE,3)
        Bbar(5,3:3*NNODE:3)   = dNdy(1:NNODE,1)
        Bbar(6,2:3*NNODE-1:3) = dNdy(1:NNODE,3)
        Bbar(6,3:3*NNODE:3)   = dNdy(1:NNODE,2)
        Bbar2 = 0.d0
        Bbar2(1,1:3*NNODE-2:3) = dNbardy(1:NNODE,1)-dNdy(1:NNODE,1)
        Bbar2(2,1:3*NNODE-2:3) = dNbardy(1:NNODE,1)-dNdy(1:NNODE,1)
        Bbar2(3,1:3*NNODE-2:3) = dNbardy(1:NNODE,1)-dNdy(1:NNODE,1)
        Bbar2(1,2:3*NNODE-1:3) = dNbardy(1:NNODE,2)-dNdy(1:NNODE,2)
        Bbar2(2,2:3*NNODE-1:3) = dNbardy(1:NNODE,2)-dNdy(1:NNODE,2)
        Bbar2(3,2:3*NNODE-1:3) = dNbardy(1:NNODE,2)-dNdy(1:NNODE,2)
        Bbar2(1,3:3*NNODE:3)   = dNbardy(1:NNODE,3)-dNdy(1:NNODE,3)
        Bbar2(2,3:3*NNODE:3)   = dNbardy(1:NNODE,3)-dNdy(1:NNODE,3)
        Bbar2(3,3:3*NNODE:3)   = dNbardy(1:NNODE,3)-dNdy(1:NNODE,3)
        Bbar2 = (1.d0/3.d0)*Bbar2
        Bbar = Bbar + Bbar2
        
        !Assemble Bstar
        Bstar = 0.d0
        Bstar(1,1:3*NNODE-2:3) = dNdy(1:NNODE,1)
        Bstar(2,2:3*NNODE-1:3) = dNdy(1:NNODE,2)
        Bstar(3,3:3*NNODE:3)   = dNdy(1:NNODE,3)
        Bstar(4,1:3*NNODE-2:3) = dNdy(1:NNODE,2)
        Bstar(5,2:3*NNODE-1:3) = dNdy(1:NNODE,1)
        Bstar(6,1:3*NNODE-2:3) = dNdy(1:NNODE,3)
        Bstar(7,3:3*NNODE:3)   = dNdy(1:NNODE,1)
        Bstar(8,2:3*NNODE-1:3) = dNdy(1:NNODE,3)
        Bstar(9,3:3*NNODE:3)   = dNdy(1:NNODE,2)
        Bstar2 = 0.d0
        Bstar2(1,1:3*NNODE-2:3) = dNbardy(1:NNODE,1)-dNdy(1:NNODE,1)
        Bstar2(2,1:3*NNODE-2:3) = dNbardy(1:NNODE,1)-dNdy(1:NNODE,1)
        Bstar2(3,1:3*NNODE-2:3) = dNbardy(1:NNODE,1)-dNdy(1:NNODE,1)
        Bstar2(1,2:3*NNODE-1:3) = dNbardy(1:NNODE,2)-dNdy(1:NNODE,2)
        Bstar2(2,2:3*NNODE-1:3) = dNbardy(1:NNODE,2)-dNdy(1:NNODE,2)
        Bstar2(3,2:3*NNODE-1:3) = dNbardy(1:NNODE,2)-dNdy(1:NNODE,2)
        Bstar2(1,3:3*NNODE:3)   = dNbardy(1:NNODE,3)-dNdy(1:NNODE,3)
        Bstar2(2,3:3*NNODE:3)   = dNbardy(1:NNODE,3)-dNdy(1:NNODE,3)
        Bstar2(3,3:3*NNODE:3)   = dNbardy(1:NNODE,3)-dNdy(1:NNODE,3)
        Bstar2 = (1.d0/3.d0)*Bstar2
        Bstar  = Bstar + Bstar2
        
        !Compute Sigma - Geometric Stiffness Term
        S    = 0.d0
        Svec = 0.d0
        Smat = 0.d0
        Pvec = 0.d0
        Pmat = 0.d0
        Sigma = 0.d0
        S = matmul(Kirch,transpose(dNdy(1:NNODE,1:3)))
        do i = 1,NNODE
           Pvec(1:3*NNODE) = reshape(spread(transpose(dNdy(i:i,1:3)),
     1                              dim=2,ncopies=NNODE),(/3*NNODE/))
           Pmat(3*i-2:3*i,1:3*NNODE) = spread(Pvec,dim=1,ncopies=3)
           Svec(1:3*NNODE) = reshape(spread(S(1:3,i:i),
     1                              dim=2,ncopies=NNODE),(/3*NNODE/))
           Smat(3*i-2:3*i,1:3*NNODE) = spread(Svec,dim=1,ncopies=3)
        end do
        Sigma(1:3*NNODE,1:3*NNODE) = Pmat(1:3*NNODE,1:3*NNODE)*
     1                            transpose(Smat(1:3*NNODE,1:3*NNODE))
        
        !Compute Q - Geometric Stiffness Term
        Q    = 0.d0
        Pvec = 0.d0
        Pmat = 0.d0
        do i = 1,NNODE
            Pvec(1:3*NNODE) = reshape(spread(transpose(dNdy(i:i,1:3)),
     1                                dim=2,ncopies=NNODE),(/3*NNODE/))
            Pmat(3*i-2:3*i,1:3*NNODE) = spread(Pvec,dim=1,ncopies=3)
        end do
        Q(1:3*NNODE,1:3*NNODE) = Pmat(1:3*NNODE,1:3*NNODE)*
     1                             transpose(Pmat(1:3*NNODE,1:3*NNODE))
        
        
        RHS(1:3*NNODE,1) = RHS(1:3*NNODE,1)
     1   - matmul(transpose(Bbar(1:6,1:3*NNODE)),stress(1:6))*
     2                                          w(kint)*determinant

        AMATRX(1:3*NNODE,1:3*NNODE) = AMATRX(1:3*NNODE,1:3*NNODE)
     1  + (matmul(transpose(Bbar(1:6,1:3*NNODE)),
     2     matmul(D,matmul(G,Bstar(1:9,1:3*NNODE))))
     3   - Sigma(1:3*NNODE,1:3*NNODE)
     4   + ((skk/3.d0)*(P(1:3*NNODE,1:3*NNODE)+Q(1:3*NNODE,1:3*NNODE))))
     5    *w(kint)*determinant
        

        ENERGY(2) = ENERGY(2)
     1   + 0.5D0*dot_product(stress,strain)*w(kint)*determinant       ! Store the elastic strain energy

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
     1                 - N2(1:nfacenodes)*adlmag(j,1)*norm(1:3)*w(kint)! Note determinant is already in normal
            end do
        end do
      end do

      return

      END SUBROUTINE UEL_BB