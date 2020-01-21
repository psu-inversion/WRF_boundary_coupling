#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([WRF_boundary_coupling], [1.0-dev1],
    [https://github.com/psu-inversion/WRF_boundary_coupling/issues/new],
    [WRF_boundary_coupling], [https://psu-inversion.github.io/WRF_boundary_coupling/])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIRS([m4])
AC_CONFIG_SRCDIR([module_CT_lib.f90])
AM_INIT_AUTOMAKE([-Wall -Werror foreign dist-bzip2 dist-xz dist-zip])

AX_IS_RELEASE([dash-version])
AX_CHECK_ENABLE_DEBUG()

# Checks for programs.
AC_PROG_FC

# Checks for libraries.
AX_LIB_NETCDF4()

if test "x${with_netcdf4_fortran}" = "xno"
then
    AC_MSG_ERROR([Can't find NetCDF4 Fortran bindings])
fi

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AX_FCFLAGS_WARN_ALL()

# Checks for library functions.

# Check which programs to build
if test "x${ax_is_release}" = "xno"
then
  AX_CONFIG_FEATURE_DEFAULT_ENABLED
else
  AX_CONFIG_FEATURE_DEFAULT_DISABLED
fi

AX_CONFIG_FEATURE([geos-chem], [turns on/off GEOS-Chem support],
	          BUILD_GEOSCHEM, [Define if you want GEOS-Chem support])
AX_CONFIG_FEATURE([carbontracker], [turns on/off CarbonTracker support],
                  BUILD_CARBONTRACKER, [Define if you want CarbonTracker support])

if test "x${ax_config_feature_geos_chem}" = xno -a "x${ax_config_feature_carbontracker}" = xno
then
    AC_MSG_WARN([No programs selected for build])
fi

AM_CONDITIONAL([GEOSCHEM], [test "x${ax_config_feature_geos_chem}" = xyes])
AM_CONDITIONAL([CARBONTRACKER], [test "x${ax_config_feature_carbontracker}" = xyes])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
