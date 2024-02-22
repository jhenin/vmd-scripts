# Make a bond between two atoms rotatable using the keyboard or the rotate_group procedure

proc make_bond_rotatable { a1 a2 { selText protein } {molid top} } {

  # First, split the selection into fragments separated by selected bond
  # We will rotate the smaller fragment

  # Use bond distance heuristic to follow molecular graph. We could use VMD's defined bonds instead.
  set bond_cutoff 1.6

  # Extra stopping criterion to avoid infinite loops
  set max_iter 1000

  set n1 0
  set part1 [atomselect $molid "$selText and not index $a2 and within $bond_cutoff of index $a1"]
  set n2 0
  set part2 [atomselect $molid "$selText and not index $a1 and within $bond_cutoff of index $a2"]

  # Recursively discover two halves until either one converges
  set i 0
  while { $i < $max_iter && $n1 < [$part1 num] && $n2 < [$part2 num]} {
    set n1 [$part1 num]
    set new_part1 [atomselect $molid "$selText and not index $a2 and within $bond_cutoff of index [$part1 list]"]
    $part1 delete
    set part1 $new_part1

    set n2 [$part2 num]
    set new_part2 [atomselect $molid "$selText and not index $a1 and within $bond_cutoff of index [$part2 list]"]
    $part2 delete
    set part2 $new_part2
    incr i
  }
  # For debugging: set beta field
  # [atomselect $molid all] set beta 0
  # $part1 set beta 1
  # $part2 set beta 2
  # puts "$part1 $part2"
  # $part1 global
  # $part2 global

  # Persistent atomselects for individual atoms, used by rotate_group
  set sela1 [atomselect $molid "index $a1"]
  set sela2 [atomselect $molid "index $a2"]
  $sela1 global
  $sela2 global

  # First partition to converge is smallest, make that one rotate
  if { $n1 == [$part1 num] } {
    $part1 global
    puts "Done. Use arrow keys to rotate."
    user add key Left "rotate_group $sela2 $sela1 $part1 -5"
    user add key Right "rotate_group $sela2 $sela1 $part1 5"
  } elseif { $n2 == [$part2 num] } {
    $part2 global
    puts "Done. Use arrow keys to rotate."
    user add key Left "rotate_group $sela1 $sela2 $part2 -5"
    user add key Right "rotate_group $sela1 $sela2 $part2 5"
  } else {
    puts "No convergence after $i iterations"
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