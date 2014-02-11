with Ada.Text_IO; use Ada.Text_IO;
with Alignment;
use Alignment;


procedure Calculate is
   Sec1: DNA:="ACACACTATTTT";
   Sec2: DNA:="AGCACACATGT";
   X: SW_DNA(Sec1'Length+Sec2'Length);

begin

   X:=SW_Al_DNA(Sec1,Sec2);
   New_Line;
   Print_DNA(X.Frst);

   New_Line;
   Print_DNA(X.Snd);

end Calculate;
