{ pkgs ? import <nixpkgs> {} }:

pkgs.mkShell {
  name = "cpp project";
  buildInputs = with pkgs; [
     cmake
     ninja
     pkg-config
     swig
     pcre
     bison
     flex
     openssl
     mpi
  ];
  
  shellHook = ''
   # Some bash command and export some env vars.
   export FOO=BAR
   echo "Starting new shell";
  '';
}