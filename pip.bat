echo off
where anaconda > temp.txt
set /p VAR=<temp.txt
echo Anaconda path is %VAR%
echo[
set MOD=%VAR:anaconda.exe=pip.exe%
REM echo Pip path is %MOD%
%MOD% install wxpython==4.0.7
del temp.txt
pause
