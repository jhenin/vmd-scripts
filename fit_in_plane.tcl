# Rotational fit around z axis only
# used to prevent spin/swivel motion of a membrane protein,
# without changing membrane orientation (no tilt)
#
# Jérôme Hénin <jerome.henin@cnrs.fr>

# To do before running:
# center alpha
# qwrap

set r [atomselect top "alpha and structure E H" frame 0]
set a [atomselect top "index [$r list]"] ;#same atoms as ref
set all [atomselect top all]

set nf [molinfo top get numframes]

for {set f 1} { $f < $nf} {incr f} {
    $a frame $f
    set m [measure fit $a $r]
    # Sine and cosine of z rotation are in first column of 4x4 matrix
    set m00 [lindex $m 0 0]
    set m01 [lindex $m 0 1]
    # Negative sign to reverse rotation
    # Use atan2 for safety
    set alpha [expr {-180.0 / $M_PI * atan2($m01, $m00)}]

    $all frame $f
    # Apply reverse
    $all move [transaxis z $alpha deg]
    # Check angles: second run of this script should find nearly zero
    puts "$f $alpha"
}

