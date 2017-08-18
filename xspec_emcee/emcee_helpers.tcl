# helper functions for emcee_xspec

set HSTART ">EMCEE>"
set HEND "<EMCEE<"
set EMCEE_DEBUG 0

# startup
proc emcee_startup {} {
    global EMCEE_DEBUG
    puts "Running emcee startup script"
    autosave off
    if { ! $EMCEE_DEBUG } {
        chatter 1
    }
    query no
}

proc emcee_tcloutr { args } {
    global HSTART HEND
    set res [eval tcloutr $args]
    puts "$HSTART$res$HEND"
}

# get statistic
proc emcee_statistic { } {
    global HSTART HEND
    puts "$HSTART[tcloutr stat]$END"
}

proc emcee_batch { pars } {
    global HSTART HEND

    set stats ""
    foreach ps $pars {
        foreach subp $ps {
            eval newpar $subp
        }
        lappend stats [tcloutr stat]
    }

    puts "$HSTART$stats$END"
}

# loop taking parameters and returning results
# exits when quit is entered or stdin closes
proc emcee_loop { } {
    global EMCEE_DEBUG HSTART HEND

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
	    puts "$HSTART -1 $HEND"
	    continue
        } else {
            eval $line
        }
    }
}
