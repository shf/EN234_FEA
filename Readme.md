# EN234 for Linux

This fork of EN234 is slightly modified to be compatible with WSL and can be compiled entirely with gfortran. A bash file is added to compile the codes based on their dependencies and make the executable file in Output_File folder. The results can be visualized in Tecplot.

# Running examples

* Go to the `src` directory, edit the `main.f90` file, and uncomment the input and output files you would like to run. 
* In `input_files`, you can choose relevant input models.
* In `user_codes\src`, you need to develop the user material or user element subroutine. Then, change the name of subroutine to `UEL`, `UMAT`, or `VUMAT`, based on your analysis and check whether any other subroutine has the same name.
* Build and run your code using `bash.sh`. 

    - Change root folder in `main.f90` to the root folder of your project, then:

`chmod +x bash.sh`

`.\bash.sh`

* Results will be written in `Output_Files`. The results will include an ASCII file with a `.out` extension, which contains error messages and information; and one or more files with a `.dat` extension, which can be read and plotted with Tecplot. 

## TODO

1. Make the directories correct
2. Add 2d user element to the code 