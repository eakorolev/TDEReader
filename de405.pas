{* JDEread: Java reader for JPL DE/LE ephemerides.
 * Copyright (C) 2009,2013  Peter Hristozov. All rights reserved.
 *
 * This file is part of JDEread.
 *
 * JDEread is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *}
unit DE405;

interface

uses
  Classes, SysUtils, Math, DEheader;

type
  {**
   * This class sets constants and read from ASCII file coefficients of Chebyshev
   * polynomials for JPL ephemeris DE405. Class is the successor of DEheader. All
   * necessary files can be found in :
   * "ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/de405/"
   *
   * @author Peter Hristozov E-mail: peterhri@hotmail.com
   * @version 1.3 2011 Jun 15.
   *}
  TDE405 = class(TDEheader)
    private
      {*
       * Array with starting dates for each file with chebychev coefficients.
       *}
      startFileDates: TDEExtArr;
      {*
       * Array with the names of each file with chebychev coefficients.
       *}
      fileNames: TDEStrArr;
    public
      {*
       * Constructor which sets the main parameters of the ephemeris.
       *}
      constructor Create;
      {*
       * Method to read the DE405 ASCII ephemeris file corresponding to 'jultime'.
       * The start and final dates of the ephemeris file are set, as are the
       * Chebyshev polynomials coefficients for Mercury, Venus, Earth-Moon, Mars,
       * Jupiter, Saturn, Uranus, Neptune, Pluto, Geocentric Moon, Sun, nutations
       * and lunar librations.
       *
       * @param jultime
       *            Julian date for calculation.
       * @return true if all reading procedures is OK and 'jultime' is in properly
       *         interval.
       *}
      function readEphCoeff(julTime: Extended): Boolean; override;
  end;

implementation

constructor TDE405.Create;
const
  _startFileDates: array[0..30] of Extended = ( 2305424.5, 2312752.5, 2320048.5,
    2327344.5, 2334640.5, 2341968.5, 2349264.5, 2356560.5, 2363856.5, 2371184.5,
    2378480.5, 2385776.5, 2393104.5, 2400400.5, 2407696.5, 2414992.5, 2422320.5,
    2429616.5, 2436912.5, 2444208.5, 2451536.5, 2458832.5, 2466128.5, 2473456.5,
    2480752.5, 2488048.5, 2495344.5, 2502672.5, 2509968.5, 2517264.5, 2524624.5);
  _fileNames: array[0..30] of String = ( 'ascp1600.405', 'ascp1620.405',
    'ascp1640.405', 'ascp1660.405', 'ascp1680.405', 'ascp1700.405', 'ascp1720.405',
    'ascp1740.405', 'ascp1760.405', 'ascp1780.405', 'ascp1800.405', 'ascp1820.405',
    'ascp1840.405', 'ascp1860.405', 'ascp1880.405', 'ascp1900.405', 'ascp1920.405',
    'ascp1940.405', 'ascp1960.405', 'ascp1980.405', 'ascp2000.405', 'ascp2020.405',
    'ascp2040.405', 'ascp2060.405', 'ascp2080.405', 'ascp2100.405', 'ascp2120.405',
    'ascp2140.405', 'ascp2160.405', 'ascp2180.405', '' );
begin
  startFileDates := _startFileDates;
  fileNames := _fileNames;

  inherited Create;

  {*
   * Numbers per interval.<br>
   * NCOEFF-2. NCOEFF is from file header.405.
   *}
  fNumbersPerInterval := 1016;

  { GROUP 1010 from file header.405 }
  fDEnomber := 405;

  { GROUP 1030 from file header.405 }
  {*
   * Start epoch for all ephemeris in Julian Days.
   *}
  fStartEpoch := 2305424.50;
  {*
   * End epoch for all ephemeris in Julian Days.
   *}
  fFinalEpoch := 2524624.50;
  {*
   * Ephemerides files are broken into intervals of length
   * "interval duration", in [days].
   *}
  fIntervalDuration := 32;

  { GROUP 1040 - 1041 from file header.405 }
  {*
   * Speed of light in [km/sec].
   *}
  fCLight := 0.299792457999999984e+06;
  {*
   * Length of an A.U., in [km].
   *}
  fAU := 0.149597870691000015e+09;
  {*
   * Earth - Moon mass ratio.
   *}
  fEMrat := 0.813005600000000044e+02;
  {*
   * Mass parameter for Mercury GM in [AU^3/day^2].
   *}
  fGM1 := 0.491254745145081187e-10;
  {*
   * Mass parameter for Venus GM in [AU^3/day^2].
   *}
  fGM2 := 0.724345248616270270e-09;
  {*
   * Mass parameter for Earth-Moon barycenter GM in [AU^3/day^2].
   *}
  fGMB := 0.899701134671249882e-09;
  {*
   * Mass parameter for Mars GM in [AU^3/day^2].
   *}
  fGM4 := 0.954953510577925806e-10;
  {*
   * Mass parameter for Jupiter GM in [AU^3/day^2].
   *}
  fGM5 := 0.282534590952422643e-06;
  {*
   * Mass parameter for Saturn GM in [AU^3/day^2].
   *}
  fGM6 := 0.845971518568065874e-07;
  {*
   * Mass parameter for Uranus GM in [AU^3/day^2].
   *}
  fGM7 := 0.129202491678196939e-07;
  {*
   * Mass parameter for Neptune GM in [AU^3/day^2].
   *}
  fGM8 := 0.152435890078427628e-07;
  {*
   * Mass parameter for Pluto GM in [AU^3/day^2].
   *}
  fGM9 := 0.218869976542596968e-11;
  {*
   * Mass parameter for Sun GM in [AU^3/day^2].
   *}
  fGMS := 0.295912208285591095e-03;

  {*
   * A new instance of variable for Chebyshev polynomials coefficients
   * with dimension (number of intervals * numbers_per_interval+1).
   *}
  SetLength(fEphemerisCoefficients, 233681);
end;

function TDE405.readEphCoeff(julTime: Extended): Boolean;
var
  mantissa1, mantissa2, exponent, i, records, j: Integer;
  fileName, line, num: String;
  rFile: TextFile;
begin
  Result := False;

  if (julTime < fStartEpoch) or (julTime >= fFinalEpoch) then Exit;

  if (julTime < fEphemerisDates[1]) or (julTime >= fEphemerisDates[2]) then begin
    try
      // Select the proper ephemeris file
      for i := 0 to Length(startfiledates) - 2 do begin
	if (julTime >= startFileDates[i]) and (julTime < startFileDates[i + 1]) then begin
	    fEphemerisDates[1] := startFileDates[i];
	    fEphemerisDates[2] := startFileDates[i + 1];
	    fileName := fileNames[i];
	    records := floor((fEphemerisDates[2] - fEphemerisDates[1]) / intervalDuration);
        end;
      end;
      fileName := fPathEph + fileName;

      AssignFile(rFile, fileName);
      Reset(rFile);

      // Read each record in the file
      for j := 1 to records do begin
        // read line 1 and ignore
        ReadLn(rFile, line);

	// read lines 2 through 341 and parse as appropriate
	for i := 2 to 341 do begin
          ReadLn(rFile, line);
	  if i > 2 then begin
            // parse first entry
            num := Copy(Line, 2, 25);
            num[22] := 'E';
            fephemeriscoefficients[(j - 1) * fnumbersperinterval + (3 * (i - 2) - 1)] := StrToFloat(num);
	  end;
          if (i > 2) and (i < 341) then begin
            // parse second entry
            num := Copy(Line, 28, 25);
            num[22] := 'E';
            fephemeriscoefficients[(j - 1) * fnumbersperinterval + (3 * (i - 2))] := StrToFloat(num);
          end;
          if i < 341 then begin
            // parse third entry
            num := Copy(Line, 54, 25);
            num[22] := 'E';
            fephemeriscoefficients[(j - 1) * fnumbersperinterval + (3 * (i - 2) + 1)] := StrToFloat(num);
          end;
        end;
      end;
      CloseFile(rFile);
      Result := true;
    except
      on E: Exception do
        WriteLn('Error '+E.ClassName+': '+E.Message);
    end;
  end else
    Result := True;
end;

end.

