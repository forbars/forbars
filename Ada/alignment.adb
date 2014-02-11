
package body Alignment is

   function weight_SW (One, Two: DNA_nuc) return Integer is
   begin
      if One=Two and One/='-' then
         return 2;
      elsif One=Two and One='-' then
         return 0;
      else
         return -1;
      end if;
   end weight_SW;

   function Get_Maximum (Matr : H_matrix) return Pair_indx is
      Maximum : Pair_indx;
   begin
      Maximum.Fst:=1;
      Maximum.Sd:=1;
      Maximum.Val:=Matr(1,1);
      for I in Matr'Range(1) loop
         for K in Matr'Range(2) loop
             if Matr (I,K) >= Maximum.Val then
               Maximum.Fst :=I;
               Maximum.Sd :=K;
               Maximum.Val:=Matr(I,K);
            end if;
         end loop;
     end loop;
     return Maximum;
   end Get_Maximum;

   function SW_Al_DNA(First_DNA, Second_DNA: DNA) return SW_DNA is
      X: SW_DNA(First_DNA'Length + Second_DNA'Length);
      H: H_matrix(1 .. (First_DNA'Length+1), 1 .. (Second_DNA'Length+1));
      Start: Pair_indx;
      N: Positive;


   begin

      for I in X.Frst'Range loop
         X.Frst(I):='-';
         X.Snd(I):='-';
      end loop;

      for I in Positive range 1 ..  First_DNA'Length+1 loop
         for K in Positive range 1 ..  Second_DNA'Length+1 loop

            if I=1 or K=1 then
               H(I,K):=0;
            else
               H(I,K):=Integer'Max(Integer'Max(Integer'Max(0,H(I-1,K-1)+weight_SW(First_DNA(I-1),Second_DNA(K-1))), H(I-1,K)+weight_SW(First_DNA(I-1),'-')),
                                  H(I,K-1)+weight_SW('-',Second_DNA(K-1)));
            end if;

         end loop;
      end loop;

      Start:= Get_Maximum(H);


      N:=First_DNA'Length + Second_DNA'Length;
      X.Frst(N):=First_DNA(Start.Fst-1);
      X.Snd(N):=Second_DNA(Start.Sd-1);

      while Start.Fst > 1 and Start.Sd >1 loop

         if Start.Fst>1 and Start.Sd >1 then
           N:=N-1;
           if H(Start.Fst-1,Start.Sd-1)>= H(Start.Fst-1,Start.Sd) and H(Start.Fst-1,Start.Sd-1)>= H(Start.Fst,Start.Sd-1) then
              Start.Fst:=Start.Fst-1;
              Start.Sd:= Start.Sd-1;
               X.Frst(N):=First_DNA(Start.Fst);
               X.Snd(N):=Second_DNA(Start.Sd);
           elsif H(Start.Fst-1,Start.Sd)>= H(Start.Fst-1,Start.Sd-1) and H(Start.Fst-1,Start.Sd)>= H(Start.Fst,Start.Sd-1) then
              Start.Fst:=Start.Fst-1;
              X.Frst(N):=First_DNA(Start.Fst);
              X.Snd(N):='-';
           else
              Start.Sd:=Start.Sd-1;
              X.Frst(N):='-';
              X.Snd(N):=Second_DNA(Start.Sd);
            end if;
         elsif Start.Fst>1 and Start.Sd = 1 then
            N:=N-1;
            Start.Fst:=Start.Fst-1;
            X.Frst(N):=First_DNA(Start.Fst);
         else
            N:=N-1;
            Start.Sd:=Start.Sd-1;
            X.Snd(N):=Second_DNA(Start.Sd);
         end if;
      end loop;

      return X;
   end SW_Al_DNA;

   procedure Print_DNA (X: in DNA) is
   begin
      for I in X'Range loop
         DNA_IO.Put(X(I));
      end loop;
   end Print_DNA;

end Alignment;
