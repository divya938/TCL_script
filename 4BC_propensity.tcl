set field_value 0.2
set fname ../dump.efield_Z_$field_value
#set pname ../movie/corr.psf


# check for presence of coordinate file
if {! [file exists $fname]} {
   vmdcon -error "Required file '$fname' not available. Exiting..."
   quit
}

# check for presence of connectivity file
#if {! [file exists $pname]} {
#   vmdcon -error "Required file '$pname' not available. Exiting..."
#   quit
#}

### This chunk is to emulate FORTRAN floating point format ###
proc fortranformat {fmt value} {
    set f [format ${fmt} ${value}]
    regexp {%(\d+)} ${fmt} -> maxwidth
    if {[string length ${f}]>${maxwidth}} {
        return [string repeat * ${maxwidth}]
    } else {return ${f}}
}
### End of the chunk ###

set logic [file exist Handedness_percentage_$field_value.xvg]

if {$logic == 1} {
  set percentage_outfile [open Handedness_percentage_$field_value.xvg a+]
  set Done [exec grep -w "@" -c -v Handedness_percentage_$field_value.xvg]
  puts "$Done"
  puts "ALREADY $Done Frames are calculated"
  mol new   $fname first $Done type lammpstrj autobonds yes waitfor all
  #mol addfile $pname type psf autobonds yes waitfor all
} else {
   mol new   $fname  autobonds yes last 0 waitfor all type lammpstrj
   #mol addfile $pname type psf autobonds yes waitfor all
   set percentage_outfile [open Handedness_percentage_$field_value.xvg w]
   puts $percentage_outfile "@ title \"Percentage of Handedness\"" 
   flush $percentage_outfile
   puts $percentage_outfile "@ subtitle \"Efield 0.2 V/\\cE\\C\"" 
   flush $percentage_outfile
   puts $percentage_outfile "@ xaxis label \"Time (ns)\"" 
   flush $percentage_outfile
   puts $percentage_outfile "@ yaxis label  \"Percentage of left handedness\"" 
   flush $percentage_outfile
   puts $percentage_outfile "@ TYPE xy" 
   flush $percentage_outfile
   puts $percentage_outfile "@ view 0.15, 0.15, 0.75, 0.85" 
   flush $percentage_outfile
   puts $percentage_outfile "@ legend on" 
   flush $percentage_outfile
   puts $percentage_outfile "@ legend box on" 
   flush $percentage_outfile
   puts $percentage_outfile "@ legend loctype view" 
   flush $percentage_outfile
   puts $percentage_outfile "@ legend 0.78,0.8" 
   flush $percentage_outfile
   puts $percentage_outfile "@ legend length 2" 
   flush $percentage_outfile
   puts $percentage_outfile "@ s0 legend \"Stack\\s1\\N\"" 
   flush $percentage_outfile
   puts $percentage_outfile "@ s1 legend \"Stack\\s2\\N\"" 
   flush $percentage_outfile
   puts $percentage_outfile "@ s2 legend \"Stack\\s3\\N\"" 
   flush $percentage_outfile
   puts $percentage_outfile "@ s3 legend \"Stack\\s4\\N\"" 
   flush $percentage_outfile
   puts $percentage_outfile "@ s4 legend \"Stack\\s5\\N\"" 
   flush $percentage_outfile
   puts $percentage_outfile "@ s5 legend \"Stack\\s6\\N\"" 
   flush $percentage_outfile
   puts $percentage_outfile "@ s6 legend \"Stack\\s7\\N\"" 
   flush $percentage_outfile
   puts $percentage_outfile "@ s7 legend \"Stack\\s8\\N\"" 
   flush $percentage_outfile
   puts $percentage_outfile "@ s8 legend \"Stack\\s9\\N\"" 
   flush $percentage_outfile
   puts $percentage_outfile "@ s9 legend \"Total\"" 
   flush $percentage_outfile
   set Done 0
}

set nframe [molinfo top get numframes]
set freq 5000
set timestep 0.5
set dt [expr $freq*$timestep*0.001*0.001]
set all [atomselect top all]
set Nmols [llength [lsort -unique -integer [$all get residue]]]
set N_dihed [llength [[atomselect top {type 7}] get index]]
set N_tail  4
set Natoms [llength [[atomselect top {residue 0}] get index]]
set Nmols_per_stack [expr $Nmols/9]
puts  $Nmols_per_stack 


#for {set istack 0} {$istack < 9} {incr istack} {
#    lappend outfile [open Handedness_$istack.dat a+]
#    lappend per_outfile [open Handedness_percentage_$istack.dat a+]
#}

# Dihedral Index to work with           
set D_index {{0 5 12 13} {5 0 8 9} {3 2 20 21} {2 3 16 17} }
# Dihedral Index to work with           

for {set iframe 0} {$iframe < $nframe } {incr iframe} {
   molinfo top set frame $iframe
 
   puts -nonewline "\033\[1;32m"; #GREEN
   puts "\t \t Iframe Progress: $iframe/[expr $nframe-1]"
   puts -nonewline "\033\[0m";# Reset
   set Total_left 0
   for {set istack 0} {$istack < 9} {incr istack} {
      
       puts -nonewline "\033\[1;34m"; #BLUE
       puts "\t \t \t Stack Progress:$istack/8"
       puts -nonewline "\033\[0m";# Reset

       set stack_left_D 0
    #  set stack_right_D 0
       for {set iresid 0} {$iresid < $Nmols_per_stack} {incr iresid} {
           
           set effective_residi [expr [expr $istack*$Nmols_per_stack]+$iresid]
            
           set resid_index [[atomselect top "residue $effective_residi"] get index]
           set Effective_dihedral_index {}
           set left_D 0
         # set right_D 0
           set dummy {}
           for {set Ndihedral 0} {$Ndihedral < ${N_tail}} {incr Ndihedral} {
              
              for {set d_index 0} {${d_index} < 4} {incr d_index} {
              lappend dummy  [lindex ${resid_index} [lindex ${D_index} ${Ndihedral} ${d_index} ]]
              }
              lappend Effective_dihedral_index  ${dummy}
              set dummy {}
              # d_index loop
             }
             # Ndihedral loop

           for {set Ndihedral 0} {$Ndihedral < ${N_tail}} {incr Ndihedral} {
              set D [measure dihed [lindex ${Effective_dihedral_index} ${Ndihedral}]]

              if {[expr ${D} < 180.0 && ${D} >=90.0] || [expr ${D} >=-90.0 && ${D} < 0.0]} {
                   incr left_D
                  #} else {
                  # incr right_D
                  }
              # End of the if loop
             }
             # Ndihedral loop
            #puts [lindex $outfile $istack] "[expr {$iframe*$dt}] ${iresid} ${left_D} ${right_D} "
            #flush [lindex $outfile $istack] 
             set stack_left_D [expr ${stack_left_D}+${left_D}]   
            #set stack_right_D [expr ${stack_right_D}+${right_D}]   
           }
       # resid loop 
    #puts [lindex ${per_outfile} $istack] "[expr {$iframe*$dt}] [expr ${stack_left_D}*100.0/[expr $Nmols_per_stack*$N_tail]] "
    #flush [lindex ${per_outfile} $istack] 
    #puts [lindex $outfile $istack] " "
    #flush [lindex $outfile $istack] 
     if {${istack} == 0} {
         puts -nonewline $percentage_outfile "[ fortranformat %15.7f [expr {$iframe*$dt}]] \t [fortranformat %15.7f [expr ${stack_left_D}*100.0/[expr $Nmols_per_stack*$N_tail]]]"
         flush $percentage_outfile
     } elseif { ${istack} >= 1 && ${istack} <= 8 } {
         puts -nonewline $percentage_outfile "[fortranformat %15.7f [expr ${stack_left_D}*100.0/[expr $Nmols_per_stack*$N_tail]]] \t"
         flush $percentage_outfile
     } 
   set Total_left [expr ${Total_left}+${stack_left_D}]  
   }
   puts $percentage_outfile "[fortranformat %15.7f  [expr ${Total_left}*100.0/[expr $N_tail*$Nmols]]]"
   flush $percentage_outfile
# stack loop
}  
# frame loop

foreach x $percentage_outfile {
close $x
}

exit

