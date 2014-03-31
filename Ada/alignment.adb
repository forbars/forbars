
package body Alignment
  with SPARK_Mode
is

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

   function Get_Maximum (Matr : H_matrix; default: Integer) return Pair_indx is
      Maximum : Pair_indx;
   begin
      Maximum.Fst:=1;
      Maximum.Sd:=1;
      Maximum.Val:=default;
      if Matr'Length(1)>0 and Matr'Length(2)>0 then
         for I in Matr'Range(1) loop
           for K in Matr'Range(2) loop
             if Matr (I,K) >= Maximum.Val then
               Maximum.Fst:=I;
               Maximum.Sd:=K;
               Maximum.Val:=Matr(I,K);
            end if;
         end loop;
         end loop;
      end if;
     return Maximum;
   end Get_Maximum;

   function SW_Al_DNA(First_DNA, Second_DNA: DNA) return SW_DNA
   is

      Dist1:constant DNA_Range:=First_DNA'Length;
      Dist2:constant DNA_Range:=Second_DNA'Length;
      Length_Al: constant Positive:=(Dist1 + Dist2);
      Range1: constant Positive:=Dist1+1;
      Range2: constant Positive:=Dist2+1;
      subtype Hi is H_matrix(1..Range1,1..Range2);
      H: Hi;
      X:SW_DNA(Length_Al);
      Start: Pair_indx;
      N: Positive;

   begin

      N:=Length_Al;
      for K in 2 ..  Range2 loop
      for I in 2 ..  Range1 loop
               H(I,K):=Integer'Max(Integer'Max(Integer'Max(0,H(I-1,K-1)+weight_SW(First_DNA(I-1),Second_DNA(K-1))), H(I-1,K)-1),
                                   H(I,K-1)-1);
         end loop;
      end loop;

      Start:= Get_Maximum(H,-(Length_Al));
      X.Score:=Start.Val;

      if (Start.Fst/=1 and Start.Sd/=1) then

      X.Frst(N):=First_DNA(Start.Fst-1);
      X.Snd(N):=Second_DNA(Start.Sd-1);

      while (Start.Fst > 1 and Start.Sd >1) loop
         N:=N-1;
         if Start.Fst>1 and Start.Sd >1 then
           if H(Start.Fst-1,Start.Sd-1)>= H(Start.Fst-1,Start.Sd) and H(Start.Fst-1,Start.Sd-1)>= H(Start.Fst,Start.Sd-1) then
              Start.Fst:=Start.Fst-1;
              Start.Sd:= Start.Sd-1;
              if Start.Fst>1 and Start.Sd>1 then
               X.Frst(N):=First_DNA(Start.Fst-1);
               X.Snd(N):=Second_DNA(Start.Sd-1);
              end if;
           elsif H(Start.Fst-1,Start.Sd)> H(Start.Fst-1,Start.Sd-1) and H(Start.Fst-1,Start.Sd)> H(Start.Fst,Start.Sd-1) then
              Start.Fst:=Start.Fst-1;
              X.Frst(N):=First_DNA(Start.Fst);
              X.Snd(N):='-';
           else
              Start.Sd:=Start.Sd-1;
              X.Frst(N):='-';
              X.Snd(N):=Second_DNA(Start.Sd);
               end if;
              end if;
       end loop;

       while Start.Fst>1 and Start.Sd = 1 loop
            Start.Fst:=Start.Fst-1;
               X.Frst(N):=First_DNA(Start.Fst);
               N:=N-1;
       end loop;

       while Start.Fst=1 and Start.Sd > 1 loop
            Start.Sd:=Start.Sd-1;
            X.Snd(N):=Second_DNA(Start.Sd);
            N:=N-1;
       end loop;
     end if;
         return X;
   end SW_Al_DNA;


end Alignment;
