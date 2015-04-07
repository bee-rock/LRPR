#!/bin/sh
GRAPHS="arabic-2005 sk-2005 uk-2007-05"

for g in $GRAPHS; do
  cp ../large_scale/inout-$g-85.log $g-inout-85-1.log
  cp ../large_scale/power-$g-85.log $g-power-85-1.log
done
