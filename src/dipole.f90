




     Program XSF_symm                             ! read excitonic wave function from xsf file and make it symmetrical 
      use Module_XSF
      call getarg(1,name)                         ! read name of the file
      call read_xsf                               ! read excitonic wave function from file
      call calc_dipole_moment                     ! calculate dipole moment of exciton
     end program XSF_symm















