with Ada.Text_IO; use Ada.Text_IO;

package Alignment is

   -- definition of nucleotide type: A - adenine, G - guanine, C - cytosine, T - thymine, - is a deletion of
   -- indefinite length

   type DNA_nuc is ('A', 'T', 'C', 'G','-');

   -- type of the string containing only letters for DNA nucleotides, as
   -- they are denoted in FASTA format and defined above

   type DNA is array (Positive range <>) of DNA_nuc;

   -- type Aligned requires both aligned sequences to be of the same length,
   -- which is something desired

   package DNA_IO is new Ada.Text_IO.Enumeration_IO (DNA_nuc);

   type Aligned_DNA (Length : Natural) is
      record
         Frst: DNA(1..Length);
         Snd: DNA(1..Length);
      end record;

   type Pair_indx is
      record
         Fst: Positive;
         Sd: Positive;
         Val: Integer;
      end record;

   -- a separate type for alignment done with Smith-Waterman (SW)

   type SW_DNA is new Aligned_DNA;

   -- Aligment matrix

   type H_matrix is array (Positive range <>, Positive range <>) of Integer;

   -- weighting shceme for SW

   function Get_Maximum (Matr : H_matrix) return Pair_indx;

   function weight_SW (One, Two: DNA_nuc) return Integer;

   -- SW alignment function

   function SW_Al_DNA (First_DNA, Second_DNA : DNA) return SW_DNA;

   procedure Print_DNA (X: in DNA);

end Alignment;
