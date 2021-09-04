# Debussy v2.2 (September 2021)

## Prerequisites for DEBUSSY:
* Download and install Anaconda Python 3.8 from https://www.anaconda.com/products/individual 
* Download and install Java from www.java.com (it is suggested to download the last version available on the website:
 Java version older than April 2019 may give some issues, due to modification in the Oracle Java license).
* FOR MacOS ONLY: Install XQuartz. Download the dmg from the following website: www.xquartz.org  
 IMPORTANT:  Restart your computer after installing XQuartz

## How to install MacOSX version
*  Move into the folder DEBUSSY_v2.2 to install the programs Suite, and type on the Terminal: ./install_debussy_v2.2 <br>
At some point you will be asked for root credentials, to move some library provided by us, in the /usr/local/* and /opt/local/* folders. <br>
During this step the anaconda3 and DEBUSSY_v2.2 paths will be added at the User environmental variable. <br>
A missing python library (wxpython 4.0.7) will be installed.<br>
The installation can take some minutes.  <br>
At the end of the procedure, you will have a message “DONE!!’  and “BYE BYE” on your terminal window. 
* Inside the DEBUSSY_v2.2 folder you can find a RUN_TEST_UNIX folder, containing some files to test the Debussy workflow. Type on the Terminal:
cd RUN_TEST_UNIX 
sh drun.sh 
The run-time output of the program should appear on the Terminal ending with “Debussy simulation done”.
* Check the GUI installation, by double clicking on the debussy-suite_gui launcher in DEBUSSY_v2.2 folder. 
  If you can see on your screen the GUI starting, the installation procedure ended successfully. 


## How to install WINDOWS version
* copy cygwin64 folder under C:\  
  Any other location is not handled 
* Double click on the install_v2.2.bat file (eventually, if possible, click on the option “Run as administrator” with the right mouse button): 
  this batch file will update the variable Path of the User, by adding the C:\cygwin64\bin path and it will install the wxpython (4.0.7) module that is missing in the Anaconda package. 
	-	If the automatic setting of the environmental variables fails, you can update the variable User Path by yourself (by adding C:\cygwin64\bin\, 
    as explained here: http://www.oxfordmathcenter.com/drupal7/node/13). 
	-	If the automatic wxpython (4.0.7) installation fails, you can type from cmd.exe:
							pip3 install wxpython==4.0.7
* Open the C:\cygwin64\RUN_TEST_WIN and double click on drun.bat. 
The run-time output of the program should appear on the Terminal ending with “Debussy simulation done”. 

* Double click on the debussy-suite_GUI.bat shortcut (under C:\cygwin64\DEBUSSY_v2.2\) for launching the GUI of the Suite 
  (you can also create a desktop shortcut and associate the debussy-icon.ico available in the folder).
 If you see on your screen the GUI starting, the installation procedure ended successfully. 
