#!/bin/sh

# run xspec
# edit this for you configuration

# $1 = name of system

host=$1
if [ $host == localhost ]; then
    # avoid using ssh if running ssh on local machine
    exec xspec
else
    # modify this command to start remote process
    # here ~/.bash_profile is the script which initialises heasoft
    exec ssh -x $host "source ${HOME}/.bash_profile; exec xspec"
fi
