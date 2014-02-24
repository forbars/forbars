(* This code is heavily influenced by Benjamin Pierce's Software Foundations
It contains formalization of Smith-Waterman local alignment and some other definitions useful for the future work. This code is a subject of farther refactoring during following months. Only the case of DNA/RNA alignment is considered, general/semi-final version including proteins is scheduled to the end of the March.*)


Require Import ZArith. 
Set Implicit Arguments.

(* creating type nucleotide as a separate type describing possible nucleotides, here D stands for deletion/gap, Er - ERROR, UNKNOWN.*)

Inductive nucleotide: Type:= 
| A: nucleotide
| T: nucleotide
| G: nucleotide
| C: nucleotide
| U: nucleotide 
| D: nucleotide
| Er: nucleotide.

(* boolean equality equivalence notion for the nucleotides *)

Definition beq_nuc (a b: nucleotide): bool:=
match a with
|A => match b with
      |A => true
      |T=> false
      |G => false
      |C=> false
      |U => false
      |D => false
      |Er => false
      end
|T => match b with
      |A => false
      |T=> true
      |G => false
      |C=> false
      |U => false
      |D => false
      |Er => false
      end
|G => match b with
      |A => false
      |T=> false
      |G => true
      |C=> false
      |U => false
      |D => false
      |Er => false
      end
|C => match b with
      |A => false
      |T=> false
      |G => false
      |C=> true
      |U => false
      |D => false
      |Er => false
      end
|U => match b with
      |A => false
      |T=> false
      |G => false
      |C=> false
      |U => true
      |D => false
      |Er => false
      end
|D => match b with
      |A => false
      |T=> false
      |G => false
      |C=> false
      |U => false
      |D => true
      |Er => false
      end
|Er => match b with
      |A => false
      |T=> false
      |G => false
      |C=> false
      |U => false
      |D => false
      |Er => true
      end
end.

Inductive nucl_type: Type:=
|DNA_case
|RNA_case.


(* function that provides a complimentary nucleotide for any given *)

Definition complimentary (X: nucl_type) (n: nucleotide): nucleotide:=
match X with
|DNA_case => match n with
            | A => T
            | T => A
            | G => C
            | C => G
            | D => D
            | U => Er
            | Er => Er
            end
|RNA_case => match n with
            | A => U
            | U => A
            | G => C
            | C => G
            | T => Er
            | D => D
            | Er => Er
            end
end.

Inductive list (X: Type): Type:=
| nil: list X
| cons: X -> list X -> list X.

(* some notations to simplify writing *)

Notation "x :: l":=(cons x l) (at level 60, right associativity). Notation "[ ]" := nil.
Notation "[ x ; .. ; y ]" := (cons x .. (cons y nil) ..).

(* length of the genetic sequence *)

Fixpoint length (X: Type) (l: list X): nat:=
match l with
|nil => O
|n::l' => S (length l')
end.

(* appending two genetic sequences *) 

Fixpoint app (X: Type) (l1 l2: list X): list X:=
match l1 with
|nil => l2
|h :: t => h :: (app t l2)
end. 

Notation "x ++ y" :=(app x y)(right associativity, at level 60).

(* the head of the sequence, a helper function *)

Definition hd (default: nucleotide) (l: list nucleotide): nucleotide:= 
match l with
|nil => default
|h::l => h
end.

(* the tail of the sequence *)

Definition tl (l: list nucleotide): list nucleotide:=
match l with
|nil => nil nucleotide
|h::t => t
end.

(* adding nucleotide from the end of the sequence *)

Fixpoint snoc (l:list nucleotide) (v:nucleotide): list nucleotide:=
match l with
|nil => cons v (nil nucleotide)
|h::t => h::(snoc t v)
end. 

(* reversing the sequence *)

Fixpoint rev (l:list nucleotide): list nucleotide:=
match l with
| nil => nil nucleotide
| h::l' => (snoc (rev l') h)
end.


Inductive sprod (X: Type) : Type :=
  pair : X -> X -> sprod X.

Notation "( x , y )" := (pair x y).

(* first element of the pair *)

Definition fst {X: Type} (p : sprod X) : X := 
  match p with
  | (x,y) => x
  end.

(* second element of the pair *)

Definition snd {X: Type} (p : sprod X) : X := 
  match p with
  | (x,y) => y
  end.

(* swapping the pair of natural numbers *)

Definition swap_pair {X: Type} (p : sprod X) : sprod X := 
  match p with
  | (x,y) => (y,x)
  end.

(* dictionary type to describe the aligment 'matrix' *)

Inductive dictionary (X: Type): Type:=
|empty: dictionary X 
|record : sprod X -> Z -> dictionary X -> dictionary X.

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

Definition beq_natprod (a b: sprod nat): bool:=
match a with
|(x, y) => if (beq_nat x (fst b)) then (beq_nat y (snd b)) else false
end.

(* searching in dictionary *)

Definition natprod:=sprod nat.

Inductive zoption: Type:= 
|None: zoption
|Some: Z -> zoption.

Fixpoint find (key: natprod) (d: dictionary nat): zoption:=
match d with
| empty => None
| record k v d' => if (beq_natprod key k) then Some v else (find key d')
end.

(*a pair of gene sequences*)

Definition gene:= list nucleotide.

Definition geneprod: Type:= sprod gene.

(* weighting function, some experiments with real biological data are scheduled *)

Definition w (case: nucl_type) (l m: nucleotide): Z:=
match l with
|D => match m with
      |D => 0%Z
      |A => (-1)%Z
      |C => (-1)%Z
      |G => (-1)%Z 
      |T => (-1)%Z
      |U => (-1)%Z
      |Er => 0%Z
      end
|Er => match m with
      |D => 0%Z
      |A => (-1)%Z
      |C => (-1)%Z
      |G => (-1)%Z 
      |T => (-1)%Z
      |U => (-1)%Z
      |Er => 0%Z
      end
|A => match m with
      |D => (-1)%Z
      |A => 2%Z
      |C => (-1)%Z
      |G => (-1)%Z 
      |T => (-1)%Z
      |U => (-1)%Z
      |Er => (-1)%Z
      end
|C => match m with
      |D => (-1)%Z
      |C => 2%Z
      |A => (-1)%Z 
      |G => (-1)%Z 
      |T => (-1)%Z
      |U => (-1)%Z
      |Er => (-1)%Z
      end
|G => match m with
      |D => (-1)%Z
      |C => (-1)%Z
      |G => 2%Z 
      |A => (-1)%Z 
      |T => (-1)%Z
      |U => (-1)%Z
      |Er => (-1)%Z   
      end
|T => match case with
      |DNA_case => match m with
            |D => (-1)%Z
            |A => (-1)%Z
            |C => (-1)%Z
            |G => (-1)%Z 
            |T => 2%Z
            |Er => (-1)%Z
            |U => (-1)%Z
      end
      |RNA_case => match m with
                   |_ => (-1)%Z
                   end
     end
|U => match case with
      |RNA_case => match m with
            |D => (-1)%Z
            |A => (-1)%Z
            |C => (-1)%Z
            |G => (-1)%Z 
            |T => (-1)%Z
            |Er => (-1)%Z
            |U => 2%Z
      end
      |DNA_case => match m with
                   |_ => (-1)%Z
                   end
     end
end.

(* counting alignment score for given pair of sequences with fixed w weighting *)

Fixpoint count (case: nucl_type) (l p: gene): Z:=
match (l, p) with
|(nil, nil) => 0%Z
|(nil, (h::l)) => 0%Z
|((h::l), nil) => 0%Z
|((h::l), (h'::l')) => ((w case h h') + (count case l l'))%Z
end.

(* removing deletions from the sequence *)

Fixpoint stripD (l: gene): gene:=
match l with
|nil => nil nucleotide
|h::l => if (beq_nuc h D) then (stripD l) else (h::(stripD l))
end.

(* creating dictionary-'matrix' with first column filled with zeros *)

Fixpoint loopHI (d: dictionary nat ) (m: nat): dictionary nat:=
match m with 
|O => record (O, O) 0%Z d
|S n' => record (S n', O) 0%Z (loopHI d n')
end.

(* creating dictionary-'matrix' with first row filled with zeros *)

Fixpoint loopHJ (d: dictionary nat) (m: nat): dictionary nat:=
match m with 
|O => record (O, O) 0%Z d
|S n' => record (O, S n') 0%Z (loopHJ d n')
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

Fixpoint find' (key: natprod) (d: dictionary nat): Z:=
match d with
| empty => 0%Z
| record k v d' => if (beq_natprod key k) then v else (find' key d')
end.

Definition get_next_m_n (case: nucl_type) (d: dictionary nat) (m n: nat) (c b a: Z) (l p: gene): dictionary nat:=
record (m+1, n+1) (Z.max (c+(w case (nth (m) l D) (nth (n) p D))) (Z.max (b+(w case (nth (m) l D) D)) (Z.max (a+(w case D (nth (n) p D))) 0))) d.

Definition ha:=A::C::A::C::A::C::T::A::(nil nucleotide).
Definition ho:=A::G::C::A::C::A::C::A::(nil nucleotide).

Fixpoint get_mth (case: nucl_type) (d: dictionary nat) (m n: nat) (l p: gene): dictionary nat:=
match n with
|O => record (m,O) 0%Z d
|S n' => get_next_m_n case (get_mth case d m (n') l p) (m-1) (n') (find' (m-1,n') ((get_mth case d m (n') l p))) (find' (m,n') (get_mth case d m (n') l p)) (find' (m-1,n) d) l p
end. 


Fixpoint fillo (case: nucl_type) (d: dictionary nat) (m n: nat) (l p: gene): dictionary nat:=
match m with
|O => loopHJ d n
|S m' => get_mth case (fillo case d m' n l p) m n l p
end. 

Definition fuh:= fillo DNA_case (empty nat) ((length ho)) ((length ha)) ho ha.


(* calculating the maximal score for aligment based on H(i,j) *)

Fixpoint maxd (d: dictionary nat) (default: Z): Z:=
match d with
|empty => default
|record i v d' => (Z.max v (maxd d' default))
end.

Definition sw (case: nucl_type) (l p: gene): Z:= maxd (fillo case (empty nat) (length (stripD l)) (length (stripD p)) l p) (-1*(Z.of_nat (length l))).
 
Definition aligned (case: nucl_type) (l p: gene): Prop:= ((length l = length p)) /\ ((count case l p) = (sw case l p)).


Fixpoint maxd_ind (d: dictionary nat): natprod:=
match d with
|empty => (0,0)
|record i v d' => if (Z.eqb v (maxd d (0)%Z)) then i else maxd_ind d'
end.

Compute maxd_ind fuh.


Definition next_step (d: dictionary nat) (pi: natprod) (l p: gene): geneprod:=
match pi with
|(O,O) => (D::(nil nucleotide),D::(nil nucleotide))
|(O, S n') => (D::(nil nucleotide), (nth n' p D)::(nil nucleotide))
|(S m', O) => ((nth m' l D)::(nil nucleotide), D::(nil nucleotide))
|(S m', S n') => if (Z.eqb (Z.max (Z.max (find' (m', n') d) (find' (m', S n') d)) (find' ( S m',n') d))  (find' (m', n') d)) then ((nth m' l D)::(nil nucleotide),(nth n' p D)::(nil nucleotide)) else 
                 if (Z.eqb (Z.max (Z.max (find' (m', n') d) (find' (m', S n') d)) (find' ( S m',n') d))  (find' (S m', n') d)) then (D::(nil nucleotide),(nth m' l D)::(nil nucleotide)) else ((nth n' p D)::(nil nucleotide),D::(nil nucleotide))
end.

Definition next_step_id (d: dictionary nat) (pi: natprod) (l p: gene): natprod:=
match pi with
|(O,O) => (0,0)
|(O, S n') => (0, n')
|(S m', O) => (m', 0)
|(S m', S n') => if (Z.eqb (Z.max (Z.max (find' (m', n') d) (find' (m', S n') d)) (find' ( S m',n') d))  (find' (m', n') d)) then (m',n') else 
                 if (Z.eqb (Z.max (Z.max (find' (m', n') d) (find' (m', S n') d)) (find' ( S m',n') d))  (find' (S m', n') d)) then (S m', n') else (m', S n')
end.

Fixpoint sw_combine (d: dictionary nat) (ha: natprod) (i: nat) (l p: gene): geneprod:= 
match i with
|O => (nil nucleotide, nil nucleotide)
|S i' => ((fst (next_step d (fst ha, snd ha) l p)) ++ (fst (sw_combine d (next_step_id d (fst ha, snd ha) l p) i' l p) ), snd (next_step d (fst ha, snd ha) l p) ++ snd (sw_combine d (next_step_id d (fst ha, snd ha) l p) i' l p))
end. 

Definition sw_align (l p: gene):= sw_combine (fillo DNA_case (empty nat) ((length l)) ((length p)) l p) (maxd_ind (fillo DNA_case (empty nat) ((length l)) ((length p)) l p)) ((length l)+(length p)) l p.  

Compute sw_align ho ha.

Compute count DNA_case (fst (sw_align ho ha)) (snd (sw_align ho ha)).

Theorem lengthSW: forall l p: gene, (length (fst (sw_align l p)))=(length (fst (sw_align l p))).
Proof. intros l p. auto. Qed.

Theorem countSW: forall l p: gene, (count DNA_case l p) = (sw DNA_case l p).

(* The proof of the theorem or a similar one is expected in commit scheduled to 3th of March, after some rethinking (maybe) of aligment specification. A human-readable annotation/ documentation of the code will be added in the same commit *)