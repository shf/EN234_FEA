!    The file contains the following subrouines:

!          abq_inverse_LU                       - inverse of matrix with >3 order using LU method
!          abq_UEL_invert2D                     - inverse of 2D matrix
!          abq_UEL_invert3D                     - inverse of 3D matrix

      subroutine abq_UEL_invert3d(A,A_inverse,determinant)
  
      double precision, intent(in) :: A(3,3)
      double precision, intent(out) :: A_inverse(3,3)
      double precision, intent(out) :: determinant
  
      double precision COFACTOR(3,3)
  
!   Compute inverse and determinant of 3x3 matrix
  
      determinant =   A(1,1)*A(2,2)*A(3,3)
     1   - A(1,1)*A(2,3)*A(3,2)
     2   - A(1,2)*A(2,1)*A(3,3)
     3   + A(1,2)*A(2,3)*A(3,1)
     4   + A(1,3)*A(2,1)*A(3,2)
     5   - A(1,3)*A(2,2)*A(3,1)
  
      IF (determinant==0.d0) THEN
        write(6,*) ' Error in subroutine abq_UEL_inver3d'
        write(6,*) ' A 3x3 matrix has a zero determinant'
        stop
      endif
      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))
  
      A_inverse = transpose(COFACTOR) / determinant
  
      end subroutine abq_UEL_invert3d


      subroutine abq_UEL_invert2d(A,A_inv,det_A)
        !
        ! Returns A_inv, the inverse, and det_A, the determinant
        ! Note that the det is of the original matrix, not the
        ! inverse
        !
        implicit none
        !
        real*8 A(2,2),A_inv(2,2),det_A,det_A_inv
  
        
        
        det_A = A(1,1)*A(2,2) - A(1,2)*A(2,1)
              
        det_A_inv = 1.d0/det_A
            
        A_inv(1,1) =  det_A_inv*A(2,2)
        A_inv(1,2) = -det_A_inv*A(1,2)
        A_inv(2,1) = -det_A_inv*A(2,1)
        A_inv(2,2) =  det_A_inv*A(1,1)
  
  
        return
        end subroutine abq_UEL_invert2d

      subroutine abq_inverse_LU(Ain,A_inverse,n)  ! Compute the inverse of an arbitrary matrix by LU decomposition
          implicit none
  
          integer, intent(in)  :: n
  
          double precision, intent(in)    :: Ain(n,n)
          double precision, intent(out)   :: A_inverse(n,n)
  
          double precision :: A(n,n), L(n,n), U(n,n), b(n), d(n), x(n)
          double precision :: coeff
          integer :: i, j, k
  
          A(1:n,1:n) = Ain(1:n,1:n)
          L=0.d0
          U=0.d0
          b=0.d0
  
          do k=1, n-1
              do i=k+1,n
                  coeff=a(i,k)/a(k,k)
                  L(i,k) = coeff
                  A(i,k+1:n) = A(i,k+1:n)-coeff*A(k,k+1:n)
              end do
          end do
  
          forall (i=1:n)  L(i,i) = 1.d0
          forall (j=1:n) U(1:j,j) = A(1:j,j)
          do k=1,n
              b(k)=1.d0
              d(1) = b(1)
              do i=2,n
                  d(i)=b(i)
                  d(i) = d(i) - dot_product(L(i,1:i-1),d(1:i-1))
              end do
              x(n)=d(n)/U(n,n)
              do i = n-1,1,-1
                  x(i) = d(i)
                  x(i)=x(i)-dot_product(U(i,i+1:n),x(i+1:n))
                  x(i) = x(i)/U(i,i)
              end do
              A_inverse(1:n,k) = x(1:n)
              b(k)=0.d0
          end do
      end subroutine abq_inverse_LU