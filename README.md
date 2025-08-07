ulc-cadical contains a preprocessor for unique literal clause (ULC) and exclusive literal clause (XLC) encoding within CaDiCaL. This preprocessor is built into CaDiCaL and called after other preprocessing techniques such as lucky search, but before the CDCL loop is entered. 

The main branch of CaDiCaL is available at https://github.com/arminbiere/cadical

Build Cadical

  > cd cadical ; ./configure && make
  OR
  > sh build.sh

Run Cadical with sequential counter encoding of ULCs with alignment on alignable formulas

  > ./build/cadical <form> proof --no-binary --ulc=1 --ulctype=1 --ulcaligntype=1 --ulcalign=1 --ulcelim=0 -t <timeout> 

Releveant options: 

--ulc=[0,1]      : 1 to perform reencoding of ulcs in preprocessing

--ulcalign=[0,1] : 0 for no alignment

--ulcelim=[0,1]  : 1 for variable elimination, transforming sequential counter to order encoding 

--ulcexit=[0,1]  : 1 to exit solver after preprocessing

--ulctype=[1,3]  : 1 for ULC detection, 3 for full XLC detection

--ulcaligntype=[0,1,2] : 0 encode no matter alignment, 1 only encode alignable formula (skip preprocessing otherwise), 2 only encode alignable and independent formula
