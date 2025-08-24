default: dipole

FC = ifort
FFLAGS = 

obj/Module_XSF.o: Module_XSF.f90
	@mkdir -p obj
	$(FC) $(FFLAGS) -c Module_XSF.f90  -o obj/Module_XSF.o  -module obj
obj/dipole.o: dipole.f90
	$(FC) $(FFLAGS) -c dipole.f90      -o obj/dipole.o      -module obj


OBJS =    obj/Module_XSF.o     \
          obj/dipole.o


dipole: $(OBJS)
	$(FC) $(LNK) -o ~/bin/dipole $(OBJS) -module obj


clean:
	rm -rf ~/bin/dipole obj/*.o obj/*.mod

