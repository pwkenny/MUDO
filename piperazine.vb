#  Vector bindings for setting protonation states of piperazines
#  Note how the [$NtoCat] definition allows only one nitrogen
#  of the piperazine to protonate.
#
Csp3          [CX4]
SecAmine      [N;H1]([$Csp3])[$Csp3]
TerAmine      [N;H0]([$Csp3])([$Csp3])[$Csp3]
NtoCat        NCC[N+]
Prot1         [$SecAmine;!$NtoCat]
Prot2         [$TerAmine;!$NtoCat]
#
