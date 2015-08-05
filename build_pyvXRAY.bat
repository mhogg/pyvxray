@echo off
REM Script use to build pyvXRAY source and binary distributions
cd cython
cython -a cythonMods.pyx
cd..
REM @CALL abaqus python setup.py sdist upload
@CALL abaqus python setup.py sdist
@CALL c:\python26_32bit\python.exe setup.py bdist
del pyvXRAY\*.pyd
@CALL abaqus python setup.py bdist
