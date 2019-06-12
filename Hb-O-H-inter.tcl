mol new 2-mer-folded.gro type gro

set O-group [list [[atomselect top "name O O1 O2 O3 O4 O5 and residue 0"] get serial] [[atomselect top "name O O1 O2 O3 O4 O5 and residue 1"] get serial] ]
set H-group [list [[atomselect top "name H3 H4 H5 H68 H72 H73 and residue 0"] get serial] [[atomselect top "name H3 H4 H5 H68 H72 H73 and residue 1"] get serial] ] 
set output [open "Hb_list-O-H-inter.dat" w]

set count 0

for {set oxy-groupid 0} {${oxy-groupid} <= 1} {incr oxy-groupid} {
    for {set hyd-groupid 0} {${hyd-groupid} <= 1} {incr hyd-groupid} {
			if {${hyd-groupid} != ${oxy-groupid} } {
            for {set i 0} {$i < 6} {incr i} {
			   for {set j 0} {$j < 6} {incr j} {
	           incr count
               puts $output "ATOMS${count}=[lindex ${O-group} ${oxy-groupid} ${i}],[lindex ${H-group} ${hyd-groupid} ${j}]"
                   }
            }
			 }
	}
	}

 puts $output "SWITCH={RATIONAL R_0=0.02 D_0=0.22 NN=6 MM=12}"

close $output
exit
