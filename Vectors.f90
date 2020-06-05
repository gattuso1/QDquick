module Vectors

use Constants_au

contains

!get vector of module d and random orientation
function vector(d)

      implicit none
      real(kind=dp), dimension(3) :: vector
      real(kind=dp) :: d, phi, teta, rand1, rand2

      call random_number(rand1)
      call random_number(rand2)
      teta = rand1*2.e0_dp*PI
      phi = acos(2*rand2-1.e0_dp)

      vector(1)=0.e0_dp
      vector(2)=0.e0_dp
      vector(3)=0.e0_dp

      vector(1)= d*sqrt(1.e0_dp-cos(phi)**2)*cos(teta)

      vector(2)= d*sqrt(1.e0_dp-cos(phi)**2)*sin(teta)

      vector(3)= d*cos(phi)

!write(6,*) vector(1), vector(2), vector(3)

end function vector

!get vector of random module and random orientation inside volume ddot
function vectorin(ddot)

      implicit none
      real(kind=dp), dimension(3) :: vectorin
      real(kind=dp) :: ddot, phi, teta, rand1, rand2, rand3, dval

      call random_number(rand1)
      call random_number(rand2)
      call random_number(rand3)
      !teta = rand1*2*PI
      !phi = acos(2*rand2-1)
      rand3 = rand3*ddot

      !vectorin(1)=0
      !vectorin(2)=0
      !vectorin(3)=0

      !vectorin(1)= dval*sqrt(1-cos(phi)**2)*cos(teta)

      !vectorin(2)= dval*sqrt(1-cos(phi)**2)*sin(teta)

      !vectorin(3)= dval*cos(phi)

    teta = rand1 * 2._dp * pi
    phi = dacos(2._dp * rand2 - 1._dp)
    dval = rand3**(1._dp/3._dp)
    vectorin(1)= dval * sin(phi) * cos(teta)
    vectorin(2)= dval * sin(phi) * sin(teta)
    vectorin(3)= dval * cos(phi)

end function vectorin

!get vector of random module and random orientation outside volume ddot
function vectorout(ddot, dmax)

      implicit none
      double precision, external:: s13adf
      real(kind=8), dimension(3) :: vectorout
      real(kind=8) :: ddot,  phi, teta, rand1, rand2, dmax, rand3, PI, dval

      call random_number(rand1)
      call random_number(rand2)
      call random_number(rand3)
      teta = rand1*2*PI
      phi = acos(2*rand2-1)
      dval = ddot + rand3 * (dmax-ddot)

      vectorout(1)=0
      vectorout(2)=0
      vectorout(3)=0

      vectorout(1)= dval*sqrt(1-cos(phi)**2)*cos(teta)

      vectorout(2)= dval*sqrt(1-cos(phi)**2)*sin(teta)

      vectorout(3)= dval*cos(phi)

end function vectorout

end module Vectors
