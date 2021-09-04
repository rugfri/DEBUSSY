@echo off
rem Set JMOL_HOME to the Jmol installation directory.
rem
if "%JMOL_HOME%x"=="x" set JMOL_HOME="C:\cygwin64\DEBUSSY_v2.2\jmol-13.0.8"
java -Xmx512m -jar "%JMOL_HOME%\Jmol.jar" %1 %2 %3 %4 %5 %6 %7 %8 %9
