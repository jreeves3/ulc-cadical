ulc-cadical contains a preprocessor for unique literal clause (ULC) and exclusive literal clause (XLC) encoding within CaDiCaL. This preprocessor is built into CaDiCaL and called after other preprocessing techniques such as lucky search, but before the CDCL loop is entered. 

The main branch of CaDiCaL is available at https://github.com/arminbiere/cadical

Build Cadical

  > cd cadical ; ./configure && make
  OR
  > sh build.sh

Run Cadical with sequential counter encoding of ULCs with alignment on alignable formulas

  > ./build/cadical <form> proof --no-binary --orderencode=1 --orderencodetype=1 --orderencodealigntype=1 --orderencodealign=1 --orderencodeelim=0 -t <timeout> 

Releveant options: 

--orderencode=[0,1]      : 1 to perform order encoding
--orderencodealign=[0,1] : 0 for no alignment
--orderencodeelim=[0,1]  : 1 for variable elimination, transforming sequential counter to order encoding 
--orderencodeexit=[0,1]  : 1 to exit solver after preprocessing
--orderencodetype=[1,3]  : 1 for ULC detection, 3 for full XLC detection
--orderencodealigntype=[0,1,2] : 0 encode no matter alignment, 1 only encode alignable formula (skip preprocessing otherwise), 2 only encode alignable and independent formula
