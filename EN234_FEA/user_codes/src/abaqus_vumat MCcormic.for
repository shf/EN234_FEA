!
!    ABAQUS format user material subroutine for explicit dynamics
!
!

      subroutine vumat(nblock, ndir, nshr, nstatev, nfieldv, nprops,
     1  lanneal, stepTime, totalTime, dt, cmname, coordMp, charLength,
     2  props, density, strainInc, relSpinInc,
     3  tempOld, stretchOld, defgradOld, fieldOld,
     4  stressOld, stateOld, enerInternOld, enerInelasOld,
     5  tempNew, stretchNew, defgradNew, fieldNew,
     6  stressNew, stateNew, enerInternNew, enerInelasNew )
!
!      include 'vaba_param.inc'
!
      implicit double precision (a-h,o-z)

      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     9  defgradNew(nblock,ndir+nshr+nshr),
     1  fieldNew(nblock,nfieldv),
     2  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     3  enerInternNew(nblock), enerInelasNew(nblock)

       character*80 cmname
!
!      Local variables
!
       integer iblock

       integer :: k, ntens, nit, maxit

       double precision :: dedev(ndir+nshr)                    ! delta(eij)
       double precision :: sdevstar(ndir+nshr)                 ! Sij^n
       double precision :: sestar                              ! predictor for effective stress
       double precision :: skkstar                             ! volumetric part of stress
       double precision :: E,xnu,Y,e0,m,edot0,alp,td,H,S,Om    ! Material properties
       double precision :: f, dfde                             ! nonlinear equation and tangent
       double precision :: eplas, deplas 
       double precision :: ta, dta, C, sig0, c1, ds0, dCde     ! variables in the nonlinear equation
!
!      Conventions for storing tensors:
!
!      Deformation gradients are provided as components in the global basis
!      (ABAQUS also allows the user to define a fixed local coord system defined with *ORIENTATION)
!
!      Stresses and stretches are stored as components in a basis that 'rotates with the material'
!      These are defined as follows:
!          Let e_i be the initial basis vectors (the global ijk directions, or user-specified basis vectors)
!          Let R be the rotation tensor, defined through the polar decomposition of deformation gradient F=RU=VR
!          Then define rotated basis vectors m_i = R e_i.
!          The co-rotational components of a tensor are defined as A_ij = m_i . A . m_j
!          The components A_ij can also be interpreted as the global components of the tensor  R^T A R
!
!
!      The components of symmetric tensors (stress,strainInc and stretch) are stored as vectors in the following order
!                      For 3D problems (ndir=3,nshr=3) [11,22,33,12,23,13]
!                      For 2D problems (ndir=3,nshr=1) [11,22,33,12]
!      The components of unsymmetric tensors (defgrad and relSpinInc) are stored as vectors in the following order
!                      For 3D problems (ndir=3,nshr=3) [11,22,33,12,23,31,21,32,13]
!                      For 2D problems (ndir=3,nshr=1) [11,22,33,12,21]
!
!      The stresses are Cauchy (true) stress
!
!      nblock                   No. integration points in data block (data are provided in 'blocks' of integration points)
!                               EN234FEA always has only 1 int pt in each block
!      ndir                     No. direct tensor components
!      nshr                     No. shear tensor components
!      nstatev                  No. user defined state variables (declared in input file)
!      nprops                   No. material properties
!      lanneal                  Annealing flag - if set to 1, then stresses are zeroed, and state vars should be reset to initial state
!      stepTime                 Elapsed time in this step at end of increment
!      totalTime                Total elapsed time at end of time increment
!      dt                       time step
!      cmname                   Material name
!      coordMp(n,i)             ith coord of nth integration point
!      charLength(i)            Characteristic length of element (in EN234FEA this is (element volume)^1/3)
!      props(k)                 kth material property
!      density                  density
!      strainInc(n,i)           ith strain increment component for nth int pt.  Strain components are in a
!                               basis that rotates with material, stored in order [de11, de22, de33, de12, de23, de13]
!                               NOTE THAT SHEAR STRAIN COMPONENTS ARE NOT DOUBLED IN THE VECTOR
!      relSpinInc(n,i)          ith component of relative spin increment tensor.   The global components of relSpinInc are
!                               dW_ij - dR_ij   where dW is the skew part of displacement gradient increment, and
!                               dR_ij is the rotation increment.
!      tempOld                  Temperature at start of increment
!      stretchOld(n,i)          ith component of right stretch V in co-rotational basis (equivalent to U in global basis) at start of increment
!      defgradOld(n,i)          ith component of deformation gradient in fixed global basis at start of increment
!      fieldOld(n,i)            ith user-defined field variable at nth integration point at start of increment
!      stressOld(n,i)           ith component of Cauchy stress in co-rotational basis (equivalent to R^T sigma R in fixed basis)
!      stateOld(n,i)            ith user defined material state variable at start of increment
!      enerInternOld(n)         specific internal energy at nth integration point, start of increment
!      enerInelasOld(n)         specific irreversible dissipation at nth integration point, start of increment
!      tempNew                  Temperature at end of increment
!      stretchNew(n,i)          ith component of right stretch V in co-rotational basis (equivalent to U in global basis) at start of increment
!      defgradNew(n,i)          ith component of deformation gradient in fixed global basis at end of increment
!      fieldNew(n,i)            ith user-defined field variable at nth integration point at end of increment
!      stressNew(n,i)           ith component of Cauchy stress in co-rotational basis (equivale nt to R^T sigma R in fixed basis)
!      stateNew(n,i)            ith user defined material state variable at end of increment
!      enerInternNew(n)         specific internal energy at nth integration point, end of increment
!      enerInelasNew(n)         specific irreversible dissipation at nth integration point, end of increment
!
!       Coded for 3D problems; small stretch assumption (works approximately for finite rotations)
!

!       Material Properties
        E     = PROPS(1)
        xnu   = PROPS(2)
        Y     = PROPS(3)
        e0    = PROPS(4)
        m     = PROPS(5)
        edot0 = PROPS(6)
        S     = PROPS(7)
        H     = PROPS(8)
        td    = PROPS(9)
        Om    = PROPS(10)
        alp   = PROPS(11)

        ntens = ndir + nshr

        do k = 1,nblock ! number of blocks (each with 1 integration point)
              ! Initialize state variables
              eplas  = stateOld(k,1)
              ta     = stateOld(k,2)
              deplas = stateOld(k,3)
              
              ! 1. Compute delta(eij)
              dedev(1:ndir) = strainInc(k,1:ndir) -
     1                        sum(strainInc(k,1:ndir))/3.D0

              dedev(ndir+1:ntens) = strainInc(k, ndir+1:ntens)

              ! 2. Calculate Sij^n
              sdevstar(1:ndir) = stressOld(k,1:ndir) - 
     1                           sum(stressOld(k,1:ndir))/3.D0
              sdevstar(ndir+1:ntens) = stressOld(k,ndir+1:ntens)
              sdevstar(1:ntens) = sdevstar(1:ntens) +
     1                            E*dedev(1:ntens)/(1.D0+xnu)

              ! calcualting sestar
              sestar = dsqrt(1.5d0)*
     1          dsqrt(dot_product(sdevstar(1:ndir),sdevstar(1:ndir)) +
     2          2.d0*dot_product(sdevstar(ndir+1:ntens),
     3                           sdevstar(ndir+1:ntens)))

              ! Second part (volumetric part) of sigma_ij^(n+1)
              skkstar = sum(stressOld(k,1:ndir))
     1             + E*sum(strainInc(k,1:ndir))/(1.d0-2.d0*xnu)

              ! Elastic increment (either no stress, or no plastic strain rate)
              if (sestar/Y<1.D-09 .or. edot0==0.D0) then
                stressNew(k,1:ntens) = sdevstar(1:ntens)
                stressNew(k,1:ndir) = stressNew(k,1:ndir) + skkstar/3.d0
                stateNew(k,1) = eplas + deplas
                stateNew(k,2) = deplas
                stateNew(k,3) = ta + dta 
              endif

              err = 1.D0
              maxit = 100
              nit = 1
              tol = 1.D-06*Y
              if (deplas==0.D0) deplas = 1.D-09/Y
              write(1,*) ' Starting Iterations '
              
         do while (err>tol)

            if (deplas<1.d-88) then
              
              dta = dt - ta*(deplas/Om)**2.d0 
              ddtade = -2.d0*ta*deplas/(Om**2.d0)

            else                            
              ! Computing delta ta
              dta = (ta*(exp(-deplas/Om)-1.D0)+Om*dt/deplas)/
     1                    (1.D0+Om/deplas*exp(-deplas/Om))

              !N-R for new deplas
              !ddtade= -ta/Om
            ddtade= ((-ta/Om*exp(-deplas/Om)-Om*dt*deplas**(-2.D0))+
     1      (ta*(exp(-deplas/Om)-1.D0)+dt*Om*1.D0/deplas)*1.D0/deplas
     2      *exp(-deplas/Om))/(1.D0+Om/deplas*exp(-deplas/Om))

            endif
            !Computing C
            C = 1.D0 - exp(-((ta+dta)/td)**alp)
            dCde = exp(-((ta+dta)/td)**alp)*(alp/td)*dtade*
     1             ((ta+dta)/td)**(alpha-1.d0)

            !Computing sig0
            sig0 = Y*(1.D0+(eplas+deplas)/e0)**m

            !Last term of nonlinear equation
            c1 = sig0/S + H*C + log(deplas/(dt*edot0))

            !Nonlinear equation
            f = sestar/S - c1 - 1.5d0*E*deplas/(S*(1.d0+xnu))

          ds0 = Y*m*(1.d0+(eplas+deplas)/e0)**(m-1.d0)*1.d0/e0
          dfde=-ds0/S - H*dCde
     1                   - 1.5d0*E/(S*(1.d0+xnu)) - 1.D0/deplas

              deplas_new = deplas - f/dfde
              if (deplas_new<0.d0) then
                deplas = deplas/10.d0
              else
                deplas = deplas_new
              endif
              
              nit = nit + 1
              err = dabs(f)
         if (nit>maxit) then
           write(6,*) ' Newton iterations in VUMAT failed to converge '
           stop
         endif
         end do
         stressNew(k,1:ntens) =
     1           ( 1.d0 - 1.5d0*deplas*E/(1.d0+xnu)/sestar )
     2                           *sdevstar(1:ntens)
         stressNew(k,1:ndir) =
     1              stressNew(k,1:ndir) + skkstar/3.d0
    
         stateNew(k,1) = eplas + deplas
         stateNew(k,2) = ta + dta
         stateNew(k,3) = deplas

        end do

!
      End Subroutine VUMAT
