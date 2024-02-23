# Make a bond between two atoms rotatable using the keyboard or the rotate_group procedure

proc pick_and_rotate { } {
  global ::vmd_pick_event

  set ::rotate_a1 -1
  # remove previous trace, if any
  trace remove variable ::vmd_pick_event write pick_rotate_callback
  trace add variable ::vmd_pick_event write pick_rotate_callback

  puts "Use mouse to pick two atoms"
  mouse mode pick
}

proc pick_rotate_callback {args} {
  global vmd_pick_atom vmd_pick_mol
  set a $::vmd_pick_atom
  if { $::rotate_a1 == -1 } {
    set ::rotate_a1 $a
  } else {
    set a2 $a
    puts "Rotating around atoms: $::rotate_a1 $a2 in molecule $vmd_pick_mol"
    make_bond_rotatable $::rotate_a1 $a2 $vmd_pick_mol
    set ::rotate_a1 -1
  }
}


proc make_bond_rotatable { a1 a2 {molid top} } {

  if [info exists ::rotate_rep] {
    set r [mol repindex $molid $::rotate_rep]
    if { $r > -1 } {
      mol delrep $r $molid
    }
    unset ::rotate_rep
  }

  # First, split the selection into fragments separated by selected bond
  # We will rotate the smaller fragment

  set sel [atomselect $molid "same fragment as index $a1 $a2"]
  set bonds [$sel getbonds]

  set part1 [list $a1]; set cur_part1 $part1
  set part2 [list $a2]; set cur_part2 $part2

  # Extra stopping criterion to avoid infinite loops
  set max_iter 10000

  # Recursively discover half until either one converges
  set i 0
  while { $i < $max_iter && [llength $cur_part1] > 0 && [llength $cur_part2] > 0 } {
    lassign [extend_half $part1 $cur_part1 $a2 $bonds] part1 cur_part1
    lassign [extend_half $part2 $cur_part2 $a1 $bonds] part2 cur_part2
    incr i
  }

  set spart1 [atomselect $molid "index $part1"]
  set spart2 [atomselect $molid "index $part2"]

  # For debugging: set beta field
  # [atomselect $molid all] set beta 0
  # $spart1 set beta 1
  # $spart2 set beta 2

  # Persistent atomselects for individual atoms, used by rotate_group
  set sela1 [atomselect $molid "index $a1"]
  set sela2 [atomselect $molid "index $a2"]
  $sela1 global
  $sela2 global

  mol color ColorID 1
  mol representation Licorice 0.3 12

  # First partition to converge is smallest, make that one rotate
  if { [llength $cur_part1] == 0 } {
    $spart1 global
    puts "Molecule partitioning done. Use arrow keys to rotate."
    user add key Left "rotate_group $sela2 $sela1 $spart1 -5"
    user add key Right "rotate_group $sela2 $sela1 $spart1 5"
    mol selection "index $part1"
    mol addrep $molid
    set repid [expr [molinfo $molid get numreps] - 1]
    set ::rotate_rep [mol repname $molid $repid]

  } elseif { [llength $cur_part2] == 0 } {
    $spart2 global
    puts "Molecule partitioning done. Use arrow keys to rotate."
    user add key Left "rotate_group $sela1 $sela2 $spart2 -5"
    user add key Right "rotate_group $sela1 $sela2 $spart2 5"
    mol selection "index $part2"
    mol addrep $molid
    set repid [expr [molinfo $molid get numreps] - 1]
    set ::rotate_rep [mol repname $molid $repid]
  } else {
    puts "Molecule partitioning did not converge after $i iterations"
  }
}


# Rotate group $sel by $angle around vector from $sela1 to $sela2

proc rotate_group { sela1 sela2 sel angle } {
  set pos1 [lindex [$sela1 get {x y z}] 0]
  set pos2 [lindex [$sela2 get {x y z}] 0]
  set axis [vecnorm [vecsub $pos1 $pos2]]

  $sel moveby [vecinvert $pos2]
  $sel move [transabout $axis $angle]
  $sel moveby $pos2
}


proc extend_half { part1 cur_part1 a2 bonds } {
  set new_part1 [list]
  foreach c $cur_part1 {
    # Iterate over bonded neighbors of recent batch of added atoms
    foreach a [lindex $bonds $c] {
      if { $a != $a2 && [lsearch -sorted -integer $part1 $a] == -1 } {
        lappend new_part1 $a
      }
    }
  }
  set cur_part1 [lsort -unique -integer $new_part1]
  # Most efficient concatenation
  set part1 [lsort -unique -integer [list {*}$part1 {*}$cur_part1]]

  return [list $part1 $cur_part1]
}