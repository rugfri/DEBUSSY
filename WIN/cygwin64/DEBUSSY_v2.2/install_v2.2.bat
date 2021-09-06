echo off
setx  /M PATH "C:\cygwin64\bin\;%PATH%"
where pip3 > temp.txt
set /p VAR=<temp.txt
echo Anaconda3 path is %VAR%
echo[
REM set MOD=%VAR:pip3.exe=pip.exe%
REM echo Pip3 path is %MOD%
%VAR% install wxpython==4.0.7
del temp.txt
setx  /M PATH "C:\cygwin64\DEBUSSY_v2.2\bin\;%PATH%"
pause

