with Ada.Text_IO; use Ada.Text_IO;
with Ada.Integer_Text_IO; use Ada.Integer_Text_IO;
with Alignment;
use Alignment;

procedure Calculate is

   Sec1: DNA:="ATTAAATA";
   Sec2: DNA:="AAAAAATTTT";
   X: SW_DNA(Sec1'Length+Sec2'Length);

   package DNA_IO is new Ada.Text_IO.Enumeration_IO (DNA_nuc);
begin

   X:=SW_Al_DNA(Sec1,Sec2);
   New_Line;
   for I in X.Frst'Range loop
         DNA_IO.Put(X.Frst(I));
      end loop;

   Put(X.Score);
   New_Line;
  for I in X.Snd'Range loop
         DNA_IO.Put(X.Snd(I));
   end loop;

end Calculate;
