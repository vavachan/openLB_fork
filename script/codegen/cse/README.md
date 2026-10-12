# OpenLB Automatic Code Generation

Common subexpression elimination (CSE) is an established approach to reducing
the arithmetic complexity of computation kernels in LBM. OpenLB currently
provides CSE-optimized versions of common LBM collision operators.

## Prerequisites:

- python
  - Mako
  - sympy
- gcc
- make

(When using Nix Flake: instantiating evironment for cse via 'nix develop .#env-code-generation')

## Execute in order

[0]: create a dynamics or operator target from the desired model to optimize:
  - compile OpenLB with the basic config.mk using gcc and add one line in the end:
    FEATURES := INSPECT_DYNAMICS EXPORT_CODE_GENERATION_TARGETS
  - now recompile OpenLB, compile the example that uses the to be optimized dynamics
    or operators and run it
  - in tmp/introspectionData you will find .dynamics and .operator files
  - put the respective file with the code you want optimized in script/codegen/cse/targets/dynamics
    or script/codegen/cse/targets/operator

[1]: adjust script/codegen/cse/Makefile, line 3 "TYPE = dynamics" or "TYPE = operator"
     depending on what you want to optimize

[2]: make prepare_extraction
[3]: make extraction
[4]: make optimize_cse

[5]: find your optimized code in expressions/(dynamics or operator)/include and use it:
  - put the .cse.h file in root: src/cse/(dynamics or operator)
  - add the lines of the generated_cse.h in /include to src/cse/(dynamics or operator)/generated_cse.h

## Special remark

If your optimized code does not give the correct results because it has too many expressions
to optimize, in cse_utils.py, you can try to change
  "return block.cse(symbols=generator, optimizations=custom_opti, order='none')" to
  "return block.cse(symbols=generator)"

For testing, if you want to disable the cse code, add "#define DISABLE_CSE" to the top of your cpp file.

## Directory tree

- root: script/codegen/cse/
  - source/               # contains python source files
  - extractions/          # performing expression extraction
  - expressions/          # cse optimizing the expression tree
  - expressions/dynamics/ # contains the optimized dynamics tuples
  - expressions/operator/ # contains the optimized non-local operators
  - targets/dynamics      # put your .dynamics files here before execution
  - targets/operator      # put your .operator files here before execution
