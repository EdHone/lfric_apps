!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Apply adjoint Helmholtz operator stored as stencil of coefficients.

! @todo Kernel does not work correctly when run in parallel; fix is being worked on

module adj_apply_helmholtz_operator_kernel_mod

  use argument_mod,      only: arg_type,              &
                               GH_FIELD, GH_REAL,     &
                               GH_SCALAR, GH_LOGICAL, &
                               GH_READ, GH_WRITE,     &
                               STENCIL, CROSS2D,      &
                               CELL_COLUMN
  use constants_mod,     only: r_solver, i_def, l_def
  use fs_continuity_mod, only: W3
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  type, public, extends(kernel_type) :: adj_apply_helmholtz_operator_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                    &
         arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, W3),                   &
         arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3, STENCIL(CROSS2D)), &
         arg_type(GH_FIELD*9, GH_REAL,    GH_READ,  W3),                   &
         arg_type(GH_SCALAR,  GH_LOGICAL, GH_READ)                         &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: adj_apply_helmholtz_operator_code
  end type adj_apply_helmholtz_operator_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------

  public :: adj_apply_helmholtz_operator_code

contains

!> @brief Apply adjoint Helmholtz operator stored as stencil of coefficients
!> @param[in]     nlayers      Number of vertical levels to solve over
!> @param[in,out] y            Pressure field to apply adjoint operator to
!> @param[in,out] x            Application of the operator to the pressure field
!> @param[in]     smap_sizes   Stencil sizes
!> @param[in]     max_length   Maximum stencil branch length
!> @param[in]     smap         Stencil dofmap
!> @param[in]     Helm_C       Diagonal entry to Helmholtz matrix
!> @param[in]     Helm_N       North (j+1) entry to Helmholtz matrix
!> @param[in]     Helm_E       East (i+1) entry to Helmholtz matrix
!> @param[in]     Helm_S       South (j-1) entry to Helmholtz matrix
!> @param[in]     Helm_W       West (j-1) entry to Helmholtz matrix
!> @param[in]     Helm_U       Upper (k+1) entry to Helmholtz matrix
!> @param[in]     Helm_UU      2nd Upper (k+2) entry to Helmholtz matrix
!> @param[in]     Helm_D       Lower (k-1) entry to Helmholtz matrix
!> @param[in]     Helm_DD      2nd Lower (k-2) entry to Helmholtz matrix
!> @param[in]     limited_area Switch to use code that can handle stencils at
!!                             edges
!> @param[in]     ndf          Number of dofs per cell for all fields, should
!!                             be = 1
!> @param[in]     undf         Size of all field arrays
!> @param[in]     map          Array containing the address of the first dof in
!!                             the column
subroutine adj_apply_helmholtz_operator_code(nlayers, y, x, smap_sizes,        &
                                             max_length, smap, Helm_C,         &
                                             Helm_N, Helm_E, Helm_S, Helm_W,   &
                                             Helm_U, Helm_UU, Helm_D, Helm_DD, &
                                             limited_area, ndf, undf, map)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers, ndf, undf
  real(kind=r_solver), intent(inout) :: y(undf), x(undf)
  integer(kind=i_def), intent(in)    :: smap_sizes(4), max_length, &
                                        smap(ndf,max_length,4)
  real(kind=r_solver), dimension(undf), intent(in) :: Helm_C,          &
                                                      Helm_N, Helm_E,  &
                                                      Helm_S, Helm_W,  &
                                                      Helm_U, Helm_UU, &
                                                      Helm_D, Helm_DD
  logical(kind=l_def), intent(in) :: limited_area
  integer(kind=i_def), intent(in) :: map(ndf)

  ! Internal variables
  integer(kind=i_def) :: k

  ! Note: global only (limited_area = .false.)
  ! This is checked for in the calling code

  do k = 2, nlayers - 1
    x(map(1) + k - 1) = x(map(1) + k - 1) + Helm_D(map(1) + k) * y(map(1) + k)
    x(map(1) + k - 2) = x(map(1) + k - 2) + Helm_DD(map(1) + k) * y(map(1) + k)
  end do

  k = 1
  x(map(1) + k - 1) = x(map(1) + k - 1) + Helm_D(map(1) + k) * y(map(1) + k)

  k = nlayers - 2
  x(map(1) + k + 1) = x(map(1) + k + 1) + Helm_U(map(1) + k) * y(map(1) + k)

  do k = 0, nlayers - 3
    x(map(1) + k + 1) = x(map(1) + k + 1) + Helm_U(map(1) + k) * y(map(1) + k)
    x(map(1) + k + 2) = x(map(1) + k + 2) + Helm_UU(map(1) + k) * y(map(1) + k)
  end do

  do k = 0, nlayers - 1
    x(smap(1,1,1)+k) = x(smap(1,1,1)+k) + Helm_C(map(1)+k)*y(map(1)+k)
    x(smap(1,2,1)+k) = x(smap(1,2,1)+k) + Helm_W(map(1)+k)*y(map(1)+k)
    x(smap(1,2,2)+k) = x(smap(1,2,2)+k) + Helm_S(map(1)+k)*y(map(1)+k)
    x(smap(1,2,3)+k) = x(smap(1,2,3)+k) + Helm_E(map(1)+k)*y(map(1)+k)
    x(smap(1,2,4)+k) = x(smap(1,2,4)+k) + Helm_N(map(1)+k)*y(map(1)+k)
  end do

end subroutine adj_apply_helmholtz_operator_code

end module adj_apply_helmholtz_operator_kernel_mod
