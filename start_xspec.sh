#!/bin/sh

# run xspec
# edit this for you configuration

# $1 = name of system

host=$1
if [ $host == localhost ]; then
    exec xspec
else
    exec ssh -x $host "source /etc/profile; do_xray_profile; exec xspec"
fi
