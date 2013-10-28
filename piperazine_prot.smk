#  These SMIRKS direct amine protonation and mono-protonation of piperazines  
#  Note the ordering of the SMIRKS which ensures that secondary amine nitrogen 
#  protonates in preference to tertiary amide.  These SMIRKS must be used with 
#  the piperazine.vb  vector bindings
#
[$Prot1:1]>>[$Prot1;h2;+:1]
[$Prot2:1]>>[$Prot2;h1;+:1]
#