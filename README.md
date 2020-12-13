# AmbulanceDeployement

## Setting up the Julia Project Locally
### This was written by Will.  Let him know if anything's confusing or if you have suggestions to make this more helpful. The Gurobi section was written by Zander.
### If you're not already familiar with the Julia shell, skim through [this link][1].  Almost all of this setup takes place in the shell.  This should get you familiar with the basic commands and interface.
1. Pull this github repo onto your local machine.  Doesn't really matter where, this won't be the active package, just need the files.
2. Inside your Julia shell, create a new package called 'AmbulanceDeployment'
   * Follow [this tutorial][2] for more details on generating a package.
   * Make sure you also change the name of the package directory to 'AmbulanceDeployment.jl' (apparently it's good convention).
   * From here, you can replace the contents of the newly generated package with the contents of this repo
   * Specifically, paste in the contents of 'AmbulanceDeployment.jl-legacy' directory from this repo, replacing the original contents of the generated package.
3. Download and set up Gurobi
    * The following steps should walk you through the process, but refer to this [Gurobi Guide][5] for additional information.
    * You'll need a license.  You can get a free academic license and download the optimizer from [this link][3]
    * You'll see three different download options on the 'Download Gurobi Optimizer' page.  Select 'Gurobi Optimizer'.
    * On Windows
        * Run the installer
        * It will make you restart your computer before the changes take effect, but this should add the Gurobi environment variable to my machine.
        * When your computer is on again, you need to activate it with a license. In Windows command line, run the following command:
            * grbgetkey <License_number>
        * It will then ask you where to store the Gurobi license file. There are two recommended options:
            * Go with the default
            * Store it in your gurobi811 folder
                * This will be considered as a non default location, but will help prevent you from losing it. Also a good idea if you have multiple versions of Gurobi.
                * Add the following environment variable:
                    * GRB_LICENSE_FILE = C:\gurobi811\gurobi.lic (or wherever your gurobi folder is)
        * Make sure all environment variables all properly set:
            * GUROBI_HOME = C:\gurobi811\win64
            * PATH has C:\gurobi811\win64\bin
            * (Optional for safety) GRB_LICENSE_FILE = C:\gurobi811\gurobi.lic (or wherever your license is)
            * Note the difference between GUROBI_HOME and PATH
    * If you have a different version of Gurobi then you need to get Gurobi v8.1.1
        * Follow all the previous steps for downloading Gurobi 8.1.1 and setting environment variables
        * In Julia terminal, check your packages and versions with either one of these commands
            * sort(collect(Pkg.installed()))
            * Pkg.status()
        * If Julia is still using a different version of Gurobi, then run this command in Julia terminal
            * Pkg.add(name=”Gurobi”, version=”0.8.1”)
            * Note that quotations might not work if you copy-paste. Instead, type the quotations yourself in Julia terminal
    * Verify that Gurobi is installed and works properly
        * In Julia terminal, check your packages and versions with either one of these commands
            * sort(collect(Pkg.installed()))
            * Pkg.status()
        * Julia should be using Gurobi 8.1.1. If not, then follow the instructions above for changing versions
        * In Julia terminal, run the following commands
            * Import Pkg
            * Pkg.add(“Gurobi”) or Pkg.add(name=”Gurobi”, version=”0.8.1”)
            * Pkg.build(“Gurobi”)
            * using Gurobi
            * Gurobi
                * This should print “Gurobi” back to the terminal if correct. Otherwise, it will say it is undefined
            * Gurobi.Model
                * This should print “Gurobi.Model” back to the terminal if correct. Otherwise, it will say it is undefined
        * If all these commands run properly, then Gurobi is set up correctly and you are good to go
4. Back in the 'AmbulanceDeployment.jl' directory in your Julia shell, activate the package with the pkg 'activate .' command.  
    * Make sure to include the dot, it points activate towards the current directory's Project.toml file.  
    * Rather than calling 'import Pkg', I just typed the ']' char, which changes the shell mode to Pkg.  Personal preference, but I thought it was a little cleaner than typing 'Pkg.' every time you need to acces Pkg. 
5. Instantiate the package, calling 'Pkg.instantiate()'
    * Since we have a Manifest.toml file, this command will download all the packages declared in that manifest.
    * If all is going well here, this should download pretty much everything you need.  I did this a week ago, so my steps could be slightly out of order, but this command did most of the heavy lifting for me.
6. Now, cd into the 'tests' directory, and call 'include(runtests.jl)'
    * You may get errors telling you to add a package.  Just do whatever the message says.
    * Keep adding packages and calling include again until you've downloaded all necessary packages
    * Once you've added all necessary packages, you'll probably get errors for syntax issues.  That's good, since we're working on syntax issues rn.

* This is where I left off.  If you make it here you should be in a good place.  If you didn't make it here, let me know of any issues and we can add the solutions to the ReadMe.  If you want to keep going feel, or if you discovered anything along the way, feel free to add what you learn to this doc.
* It'd be nice if we can get it set up so that the package we created is the same thing as the repo, that way we don't have to copy/paste files between the repo and local packages.  Less room for error.  If you have ideas, have at it.  An added benefit is that future users could add our package by cloning the GitHub link. [This][4] might be a good staring point.

-Will

## How to run after installation
1. In the Julia terminal, navigate to the AmbulanceDeployment.jl-legacy folder
2. Import Pkg
3. Hit "]" to go into Pkg mode. Here, type "activate ."
4. Go back to Julia mode by pressing backspace
5. Pkg.instantiate()
6. Pkg.add(name="Gurobi", version="0.8.1")
7. Pkg.build("Gurobi")
8. using Gurobi
9. cd into src
10. Run a file with include("file_name"). For example, include("Ambulance_Deployment_experiments.jl")


[1]: https://docs.julialang.org/en/v1/stdlib/REPL/
[2]: https://julialang.github.io/Pkg.jl/v1/creating-packages/
[3]: https://www.gurobi.com/academia/academic-program-and-licenses/
[4]: http://ucidatascienceinitiative.github.io/IntroToJulia/Html/GithubIntroduction
[5]: https://www.gurobi.com/wp-content/plugins/hd_documentations/content/pdf/quickstart_windows_8.1.pdf
