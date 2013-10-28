#  These SMIRKS control proton transfer from amine cation to 
#  unprotonated nitrogen of piperazine.  Use with piperazine.vb
#  vector bindings and exhaustive enumeration
#
#  Moves proton from tertiary amine to tertiary amine
[H:8][N+:1]1([C:7])[C:2][C:3][N;$TerAmine;+0:4][C:5][C:6]1>>[N+0:1]1([C:7])[C:2][C:3][N+:4]([H:8])[C:5][C:6]1
#  Moves protom form secondary amine to secondary amine
[H:8][N+:1]1([H:7])[C:2][C:3][N;$SecAmine;+0:4]([H:9])[C:5][C:6]1>>[N+0:1]1([H:7])[C:2][C:3][N+:4]([H:8])([H:9])[C:5][C:6]1
#
