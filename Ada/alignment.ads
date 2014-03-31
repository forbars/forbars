-- this package is not fully proven with SPARK,
-- however, contracts, pos- and preconditions will ensure reliability during runtime.
-- it is still under development.

package Alignment
   with SPARK_Mode
is

   -- nucleotide type: 'A' - adenine, 'G' - guanine, 'C' - cytosine, 'T' - thymine, '-' is a deletion of
   -- indefinite length

   type DNA_nuc is ('A', 'T', 'C', 'G','-');

   -- string containing only DNA nucleotides, as
   -- they are denoted in FASTA format and defined above

   subtype DNA_Range is Positive range 1..100000; -- some arbitraty limit fot the maximum length of DNA sequence
   type DNA is array (DNA_Range range <>) of DNA_nuc with Default_Component_Value=>'-';

   type DNA_Al is array (Positive range <>) of DNA_nuc with Default_Component_Value=>'-';

   -- type Aligned requires aligned sequences to be of the same length,
   -- which is something desired

   type Aligned_DNA (Length : Positive) is
      record
         Frst: DNA_Al(1..Length);
         Snd: DNA_Al(1..Length);
         Score: Integer:=0;
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

   type H_matrix is array (Positive range <>, Positive range <>) of Integer with Default_Component_Value=>0;

   function Get_Maximum (Matr : H_matrix; default: Integer) return Pair_indx;

   -- weighting scheme for SW

   function weight_SW (One, Two: DNA_nuc) return Integer
   with Contract_Cases => ((One = Two) and (One /= '-') => weight_SW'Result= 2,
                           (One = Two) and (One = '-') => weight_SW'Result= 0,
                          others => weight_SW'Result=-1);

   -- SW alignment function

   function SW_Al_DNA (First_DNA, Second_DNA : DNA) return SW_DNA;


end Alignment;
