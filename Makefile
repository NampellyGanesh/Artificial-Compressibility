all:
	@gfortran -fdefault-real-8 -fbacktrace -c globalvariables.f90
	@gfortran -fdefault-real-8 -fbacktrace -c geometric_pros.f90
	@gfortran -fdefault-real-8 -fbacktrace -c  flux.f90
	@gfortran -fdefault-real-8 -fbacktrace -fbounds-check -fcheck=all globalvariables.o geometric_pros.o flux.o main_program.f90 -o b.out
clean:
	@rm *.o  *.mod b.out
