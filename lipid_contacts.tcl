# Count the contacts of specific lipids with a given protein residue
# Jérôme Hénin <jerome.henin@cnrs.fr>


proc lipid_contacts { resid } {

  set lip "DMPC"
  set names [list]
  set num [dict create]
  set n 14
  set cutoff 5

  for {set i 1} {$i <= $n} {incr i} {
    lappend names "C2$i"
    lappend names "C3$i"
    dict set num "C2$i" [expr $i]
    dict set num "C3$i" [expr {$i + $n}]
  }

  set sel [atomselect top "resname $lip and name $names and within $cutoff of (protein and noh and resid $resid)"]

  set nf [molinfo top get numframes]

  for { set i 1 } { $i <= [expr $n*2] } {incr i} {
    set count($i) 0
  }

  for { set f 0 } { $f < $nf } { incr f } {
    if {[expr {$f % ($nf/10)}] == 0} { puts stdout "frame $f"}
    $sel frame $f
    $sel update
    foreach name [$sel get name] {
      incr count([dict get $num $name])
    }
  }

  for { set i 1 } { $i <= $n } {incr i} {
    puts "$i $count($i) $count([expr $i+$n])" 
  }
}
