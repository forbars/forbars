(* This code is heavily influenced by Benjamin Pierce's Software Foundations
It contains formalization of Smith-Waterman local alignment and some other definitions useful for the future work. This code is a subject of farther refactoring during following months. Only the case of DNA alignment is considered, general version for DNA/RNA and proteins is scheduled to the end of the February.*)

(* for some reason I prefer to work with Z_scope being open *)

Require Import ZArith.
Open Scope Z_scope. 
Set Implicit Arguments.

(* creating type nucleotide as a separate type describing possible nucleotides, here D stands for deletion, nucleotide includes only DNA case.*)

Inductive nucleotide: Type:= 
| A: nucleotide
| T: nucleotide
| G: nucleotide
| C: nucleotide 
| D: nucleotide.

(* boolean equality equivalence notion for the nucleotides *)

Definition beq_nuc (a b: nucleotide): bool:=
match a with
|A => match b with
      |A => true
      |T=> false
      |G => false
      |C=> false
      |D => false
      end
|T => match b with
      |A => false
      |T=> true
      |G => false
      |C=> false
      |D => false
      end
|G => match b with
      |A => false
      |T=> false
      |G => true
      |C=> false
      |D => false
      end
|C => match b with
      |A => false
      |T=> false
      |G => false
      |C=> true
      |D => false
      end
|D => match b with
      |A => false
      |T=> false
      |G => false
      |C=> false
      |D => true
      end
end.

(* function that provides a complimentary nucleotide for any given *)

Definition complimentary (n: nucleotide): nucleotide:=
match n with
| A => T
| T => A
| G => C
| C => G
| D => D
end.

(* gene - type of the sequence (list) of nucleotides, which models genetic code. At this moment I prefer definition of the list without polymorphism. Will be changed in future *)

Inductive gene: Type:=
| nil: gene
| cons: nucleotide -> gene-> gene.

(* some notations to simplify writing *)

Notation "x :: l":=(cons x l) (at level 60, right associativity). Notation "[ ]" := nil.
Notation "[ x ; .. ; y ]" := (cons x .. (cons y nil) ..).

(* length of the genetic sequence *)

Fixpoint length (l: gene): nat:=
match l with
|nil => O
|n::l' => S (length l')
end.

(* appending two genetic sequences *) 

Fixpoint app (l1 l2: gene): gene:=
match l1 with
|nil => l2
|h :: t => h :: (app t l2)
end. 

Notation "x ++ y" :=(app x y)(right associativity, at level 60).

(* the head of the sequence, a helper function *)

Definition hd (default: nucleotide) (l:gene): nucleotide:= 
match l with
|nil => default
|h::l => h
end.

(* the tail of the sequence *)

Definition tl (l: gene): gene:=
match l with
|nil => nil
|h::t => t
end.

(* adding nucleotide from the end of the sequence *)

Fixpoint snoc (l:gene) (v:nucleotide): gene:=
match l with
|nil => [v]
|h::t => h::(snoc t v)
end. 

(* reversing the sequence *)

Fixpoint rev (l:gene): gene:=
match l with
| nil => nil
| h::l' => (snoc (rev l') h)
end.

(* helper function *)

Inductive natoption : Type:=
|Some: Z -> natoption
|None: natoption.

(* another helper function *)

Inductive nucleosoption : Type:=
|Som: nucleotide -> nucleosoption
|Non: nucleosoption.

(* type describing pairs of the natural numbers, later polymorphism replacement will be used *)

Inductive natprod : Type :=
  pair : nat -> nat -> natprod.

Notation "( x , y )" := (pair x y).

(* first element of the pair of natural numbers *)

Definition fst (p : natprod) : nat := 
  match p with
  | (x,y) => x
  end.


(* second element of the pair of natural numbers *)

Definition snd (p : natprod) : nat := 
  match p with
  | (x,y) => y
  end.

(* swapping the pair of natural numbers *)

Definition swap_pair (p : natprod) : natprod := 
  match p with
  | (x,y) => (y,x)
  end.

(* dictionary type to describe the aligment 'matrix' *)

Inductive dictionary: Type:=
|empty: dictionary
|record : natprod -> Z -> dictionary -> dictionary.

(* insert element of the 'matrix' *)

Definition insert (key: natprod) (value: Z) (d: dictionary): dictionary := (record key value d).

(* defining equality condition for the pair of natural numbers *)

Fixpoint beq_nat (n m : nat) : bool :=
  match n with
  | O => match m with
         | O => true
         | S m' => false
         end
  | S n' => match m with
            | O => false
            | S m' => beq_nat n' m'
            end
  end.

Definition beq_natprod (a b: natprod): bool:=
match a with
|(x, y) => if (beq_nat x (fst b)) then (beq_nat y (snd b)) else false
end.

Definition x: bool:= (beq_natprod (O, S O) (O, S O)).

Theorem about_x: x = true.
Proof. reflexivity. Qed.

(* searching in dictionary *)

Fixpoint find (key: natprod) (d: dictionary): Z:=
match d with
| empty => 0%Z
| record k v d' => if (beq_natprod key k) then v else (find key d')
end.

(*a pair of gene sequences*)

Inductive geneprod: Type:= |pairl: gene -> gene -> geneprod.

(* first and second gene sequences from gene prod *)

Definition fstl (lp: geneprod): gene:=
match lp with
|pairl l m => l
end.
 
Definition sndl (lp: geneprod): gene:=
match lp with
|pairl l m => m
end.

(* weighting function *)

Definition w (l m: nucleotide): Z:=
match l with
|D => match m with
      |D => 0%Z
      |A => (-1)%Z
      |C => (-1)%Z
      |G => (-1)%Z 
      |T => (-1)%Z
      end
|A => match m with
      |D => (-1)%Z
      |A => 2%Z
      |C => (-1)%Z
      |G => (-1)%Z 
      |T => (-1)%Z
      end
|C => match m with
      |D => (-1)%Z
      |A => (-1)%Z
      |C => 2%Z
      |G => (-1)%Z 
      |T => (-1)%Z
      end
|G => match m with
      |D => (-1)%Z
      |A => (-1)%Z
      |C => (-1)%Z
      |G => 2%Z 
      |T => (-1)%Z
      end
|T => match m with
      |D => (-1)%Z
      |A => (-1)%Z
      |C => (-1)%Z
      |G => (-1)%Z 
      |T => 2%Z
      end
end.

(* counting alignment score for given pair of sequences with fixed w weighting *)

Fixpoint count (l p: gene): Z:=
match pairl l p with
|pairl nil nil => 0%Z
|pairl nil (h::l) => 0%Z
|pairl (h::l) nil => 0%Z
|pairl (h::l) (h'::l') => (w h h') + (count l l')
end.

(* removing deletions from the sequence *)

Fixpoint stripD (l: gene): gene:=
match l with
|nil => nil
|h::l => if (beq_nuc h D) then (stripD l) else (h::(stripD l))
end.

(* creating dictionary-'matrix' with first column filled with zeros *)

Fixpoint loopHI (d: dictionary) (m: nat): dictionary:=
match m with 
|O => record (O, O) 0%Z d
|S n' => record (S n', O) 0%Z (loopHI d n')
end.

(* creating dictionary-'matrix' with first row filled with zeros *)

Fixpoint loopHJ (d: dictionary) (m: nat): dictionary:=
match m with 
|O => record (O, O) 0%Z d
|S n' => record (O, S n') 0%Z (loopHI d n')
end.

(* finding nth nucleotide in the gene sequence *)

Fixpoint nth (n:nat) (l:gene) (default:nucleotide) {struct l} : nucleotide :=
    match n, l with
      | O, x :: l' => x
      | O, other => default
      | S m, [] => default
      | S m, x :: t => nth m t default
    end.

(* filling the H(i,j) matrix by using given sequences *)

Fixpoint loopi (d: dictionary) (m n: nat) (l p: gene): dictionary:=
match m with 
|_ => match n with
     |O => record (O, m) 0%Z d
     |S n' => record (S n', m) (Z.max (Z.max (Z.max (0%Z) ((find (pred m, n') d)+(w (nth m l D) (nth (S n') p D))))
                                            ((find ((pred m), (S n')) d)+(w (nth m l D) D)))
                                           (((find (m, n') d))+(w D (nth (S n') p D)))) d
     end

end.

Fixpoint fillo (d: dictionary) (m n: nat) (l p: gene): dictionary:=
match n with
|O => loopHI d n
|S n' => loopi (fillo d m n' l p) m (S n') l p
end. 

(* calculating the maximal score for aligment based on H(i,j) *)

Fixpoint maxd (d: dictionary) (default: Z): Z:=
match d with
|empty => default
|record i v d' => (Z.max v (maxd d' default))
end.

Definition sw (l p: gene): Z:= maxd (fillo empty (length (stripD l)) (length (stripD p)) l p) (-1*(Z.of_nat (length l))).
 
Inductive aligned: geneprod -> Prop:=
|aligned0: aligned (pairl nil nil)
|alighed1: forall l p: nucleotide, aligned (pairl (l::nil) (p::nil))
|aligned2: forall l p: gene, (length l = length p) /\ ((count l p) = (sw l p)) -> aligned (pairl l p).

