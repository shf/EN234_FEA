!    The file contains the following subrouines:

!          abq_facenodes_3D                       - returns list of nodes on the face of a 3D element
!          abq_facenodes_2D                       - returns list of nodes on the face of a 2D element

      subroutine abq_facenodes_3D(nelnodes,face,list,nfacenodes)

        implicit none
  
        integer, intent (in)      :: nelnodes
        integer, intent (in)      :: face
        integer, intent (out)     :: list(*)
        integer, intent (out)     :: nfacenodes
  
      !
      !        Subroutine to return list of nodes on an element face for standard 3D solid elements
      !
  
        if (nelnodes == 4) then
          nfacenodes = 3
          if   (face == 1) list(1:3) = [1,2,3]
          if (face == 2) list(1:3) = [1,4,2]
          if (face == 3) list(1:3) = [2,4,3]
          if (face == 4) list(1:3) = [3,4,1]
        else if (nelnodes ==6) then
          nfacenodes = 3
          if (face==1) list(1:3) = [1,2,3]
          if (face==2) list(1:3) = [6,5,4]
          if (face==3) list(1:4) = [1,2,5,4]
          if (face==4) list(1:4) = [2,3,6,5]
          if (face==5) list(1:4) = [4,6,3,1]
          if (face>2) nfacenodes = 4
        else if (nelnodes == 10) then
          nfacenodes = 6
          if   (face == 1) list(1:6) = [1,2,3,5,6,7]
          if (face == 2) list(1:6) = [1,4,2,8,9,5]
          if (face == 3) list(1:6) = [2,4,3,9,10,6]
          if (face == 4) list(1:6) = [3,4,1,10,8,7]
        else if (nelnodes == 8) then
          nfacenodes = 4
          if (face==1) list(1:4) = [1,2,3,4]
          if (face==2) list(1:4) = [5,8,7,6]
          if (face==3) list(1:4) = [1,5,6,2]
          if (face==4) list(1:4) = [2,6,7,3]
          if (face==5) list(1:4) = [3,7,8,4]
          if (face==6) list(1:4) = [4,8,5,1]
        else if (nelnodes ==15) then
          nfacenodes = 6
          if (face==1) list(1:6) = [1,2,3,7,8,9]
          if (face==2) list(1:6) = [6,5,4,11,10,12]
          if (face==3) list(1:8) = [1,2,5,4,7,14,10,13]
          if (face==4) list(1:8) = [2,3,6,5,8,15,11,14]
          if (face==5) list(1:8) = [4,6,3,1,12,15,9,13]
          if (face>2) nfacenodes = 8
        else  if (nelnodes == 20) then
          nfacenodes = 8
          if (face == 1) list(1:8) = [1,2,3,4,9,10,11,12]
          if (face == 2) list(1:8) = [5,8,7,6,16,15,14,13]
          if (face == 3) list(1:8) = [1,5,6,2,17,13,18,9]
          if (face == 4) list(1:8) = [2,6,7,3,18,14,19,10]
          if (face == 5) list(1:8) = [3,7,8,4,19,15,6,11]
          if (face == 6) list(1:8) = [4,8,5,1,20,16,17,12]
        endif
  
        end subroutine abq_facenodes_3D


        subroutine abq_facenodes_2D(nelnodes,face,list,nfacenodes)

            implicit none
      
            integer, intent (in)      :: nelnodes
            integer, intent (in)      :: face
            integer, intent (out)     :: list(*)
            integer, intent (out)     :: nfacenodes
          !
          !        Subroutine to return list of nodes on an element face for standard 2D solid elements
          !
            integer :: i3(3)
            integer :: i4(4)
      
            i3(1:3) = [2,3,1]
            i4(1:4) = [2,3,4,1]
      
            if (nelnodes == 3) then
              nfacenodes = 2
              list(1) = face
              list(2) = i3(face)
            else if (nelnodes == 4) then
              nfacenodes = 2
              list(1) = face
              list(2) = i4(face)
            else if (nelnodes == 6) then
              nfacenodes = 3
              list(1) = face
              list(2) = i3(face)
              list(3) = face+3
            else if (nelnodes == 8) then
              nfacenodes = 3
              list(1) = face
              list(2) = i4(face)
              list(3) = face+4
            else if (nelnodes == 9) then
              nfacenodes = 3
              if (face==1) list(1:3) = (/1,3,2/)
              if (face==2) list(1:3) = (/3,9,6/)
              if (face==3) list(1:3) = (/9,7,8/)
              if (face==4) list(1:3) = (/7,1,4/)
            endif
      
            end subroutine abq_facenodes_2d