# wp3

Run tests with julia test/test_wps.jl

The code is written to count self and corss counts using the same functions.
Because of this it will overcount the self pairs by a factor of 6 or 2 (depending on
whether all Self-Self-Self or Self-Self-Cross).
These factors will be present in the histogram but are removed when computing the actual DDD.