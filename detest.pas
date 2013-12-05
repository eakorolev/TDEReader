{* TDEReader: Peter Hristozov's JDEread (http://sourceforge.net/projects/jderead/) ported to Pascal.
 *
 * Ported by Evgeny Korolev E-mail: mail@eakorolev.ru
 *}
unit DEtest;

interface

uses
  Classes, SysUtils, Math, DEephem, DEheader, DE405;

{**
 * Run RunDETest to make tests.
 *
 * @param TestFilePath
 *        Path to 'testpo.XXX' file
 *}
procedure RunDETest(TestFilePath: string);

implementation

procedure RunDETest(TestFilePath: string);
const
  Header: string =
    'line# -- jed -- t# c# x# ---- NASA value ----- ---- calculated ----- ----- diff ------';
  freq: byte = 100;
var
  de: TDE405;
  test: TextFile;
  line: string;
  et, nasa, calc, delta, maxDelta: extended;
  nTarg, nCtr, nCoord, Counter: integer;
begin
  Counter := 0;
  WriteLn(#9 + 'DE test started');
  if FileExists(TestFilePath) then
    try
      AssignFile(test, TestFilePath);
      Reset(test);
      WriteLn('File ' + ExtractFileName(TestFilePath) + ' is open successfully. Ready! Set! Go!');

      ReadLn(test, line);
      case Copy(line, 1, 14) of
        'DE-0405LE-0405':
        begin
          de := TDE405.Create;
          de.pathEph := ExtractFileDir(TestFilePath) + DirectorySeparator;
        end;
        else
        begin
          WriteLn('Ephemeris "' + Copy(line, 1, 14) + '" are not supported. Exit');
          Exit;
        end;
      end;

      WriteLn('Ephemeris init: ' + de.toString);

      for ncoord := 1 to 4 do
        ReadLN(test, line);
      WriteLn(header);

      ReadLn(test, line);
      if Copy(line, 1, 3) <> 'EOT' then
      begin
        WriteLN('No "EOT" header line found. Break.');
        Exit;
      end;
      WriteLn(line);

      while not EOF(test) do
      begin
        Counter := Counter + 1;
        ReadLn(test, line);

        et := StrToFloat(Copy(line, 16, 10));
        ntarg := StrToInt(Copy(line, 27, 2));
        nctr := StrToInt(Copy(line, 30, 2));
        ncoord := StrToInt(Copy(line, 33, 2));
        nasa := StrToFloat(Copy(line, 36, 21));

        calc := TDEephem.PlEph(de, et, ntarg, nctr)[ncoord];
        delta := abs(nasa - calc);

        if (ntarg = 15) and (ncoord = 3) then
          delta := delta / (0.23 * (et - 2451545.0));
        if (ntarg = 15) and (ncoord = 6) then
          delta := delta * 0.01 / (1.0 + (et - 2451545.0) / 36525.0);

        if delta > maxDelta then
          maxDelta := delta;

        if (delta > 1e-13) or IsNan(delta) then
          WriteLn(ErrOutput, Counter:5, et: 10: 1, ntarg: 3, nctr: 3, ncoord: 3,
            nasa: 22: 13, calc: 22: 13, delta: 18);

        if Counter mod 100 = 0 then
        begin
          WriteLn(Counter:5, et: 10: 1, ntarg: 3, nctr: 3, ncoord: 3,
            nasa: 22: 13, calc: 22: 13, delta: 18);
        end;
      end;
      WriteLn('Max diff: ', maxDelta);
    except
      on E: Exception do
        WriteLn('Error ' + E.ClassName + ': ' + E.Message);
    end
  else
  begin
    WriteLn('No test file found. Exit');
  end;
end;

end.
