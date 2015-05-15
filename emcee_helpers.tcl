# helper functions for emcee_xspec
# Jeremy Sanders 2012

set HDR "@@EMCEE@@"
set EMCEE_DEBUG 0

# startup
proc emcee_startup {} {
    puts "Running emcee startup script"
    autosave off
    chatter 1
}

# get list of components and parameters
proc emcee_interrogate_params {} {
    global HDR
    puts $HDR

    query no
    fit 100

    # return component information
    #  number of components * number of datagroups
    puts [ expr [tcloutr modcomp]*[tcloutr datagrp] ]

    #  loop over data groups
    for {set dg 1} {$dg <= [tcloutr datagrp]} {incr dg} {
	#  loop over components
	for {set c 1} {$c <= [tcloutr modcomp]} {incr c} {
	   puts [tcloutr compinfo $c $dg]
	}
    }

    # return parameter information
    puts [tcloutr modpar]
    for {set i 1} {$i <= [tcloutr modpar]} {incr i} {
	puts [tcloutr pinfo $i]
	puts [tcloutr plink $i]
	if { [llength [tcloutr param $i]] > 1 } {
	    puts [tcloutr param $i]
	    puts [tcloutr sigma $i]
	} else {
	    # switch parameters only return a single value and return
	    # an error for tcloutr sigma, so work around this
	    puts "[tcloutr param $i] -1 -1e30 -1e30 1e30 1e30"
	    puts "-1"
	}
    }
    puts $HDR
}

# get statistic
proc emcee_statistic { } {
    global HDR
    puts "$HDR [tcloutr stat] $HDR"
}

proc emcee_batch { pars } {
    global HDR
    set stats ""
    foreach p $pars {
        eval newpar $p
        lappend stats [tcloutr stat]
    }
    puts "$HDR $stats $HDR"
}

# loop taking parameters and returning results
# exits when quit is entered or stdin closes
proc emcee_loop { } {
    global EMCEE_DEBUG
    global HDR

    fconfigure stdin -buffering line
    fconfigure stdout -buffering line

    while { 1 } {
	set line [gets stdin]

        if { $EMCEE_DEBUG } {
	    puts "Debug: $line"
        }

        if { [string range $line 0 4] == "batch" } {
            emcee_batch [string range $line 6 end]
	} elseif { [eof stdin] } {
	    tclexit
	} elseif { $line == "quit" } {
	    tclexit
	} elseif { $line == "returnerror" } {
	    # this is evil - asked to return an error status because
	    # the parameters were originally out
	    puts "$HDR -1 $HDR"
	    continue
        } else {
            eval newpar $line
            emcee_statistic
        }
    }
}

emcee_startup
