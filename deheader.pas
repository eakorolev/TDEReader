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
unit DEheader;

interface

uses
  Classes, SysUtils;

type
  TDEExtArr = array of Extended;
  TDEExtMatr = array of array of Extended;
  TDEIntArr = array of Integer;
  TDEStrArr = array of String;
  {*
   * Base class containing constants and parameters for any JPL ephemerides and
   * the array of coefficients for Chebyshev polynomials. For each JPL ephemeris
   * is created subclass of that class, where constants and methods are redefined.
   * Current values are for JPL ephemeris DE405.
   *
   * @author Peter Hristozov E-mail: peterhri@hotmail.com
   * @version 1.3 2011 Jun 15.
   *}
  TDEheader = class
    protected
      {*
       * Numbers per interval.<br>
       * Each interval contains an interval number, length, start and end
       * jultimes, and Chebyshev coefficients. We keep only the coefficients
       * NCOEFF-2. NCOEFF is from file header.xxx.
       *}
      fNumbersPerInterval: Integer;

      { GROUP 1010 from file header.xxx }
      {*
       * Three-digit number of planetary ephemeris version.
       *}
      fDEnomber: Integer;

      { GROUP 1030 from file header.xxx }
      {*
       * Start epoch for all ephemeris in Julian Days.
       *}
      fStartEpoch: Extended;
      {*
       * End epoch for all ephemeris in Julian Days.
       *}
      fFinalEpoch: Extended;
      {*
       * Ephemerides files are broken into intervals of length
       * "interval duration", in [days].
       *}
      fIntervalDuration: integer;

      { GROUP 1040 - 1041 from file header.xxx }
      {*
       * Speed of light in [km/sec].
       *}
      fCLight: extended;
      {*
       * Length of an A.U., in [km].
       *}
      fAU: extended;
      {*
       * Earth - Moon mass ratio.
       *}
      fEMrat: extended;
      {*
       * Mass parameter for Mercury GM in [AU^3/day^2].
       *}
      fGM1: extended;
      {*
       * Mass parameter for Venus GM in [AU^3/day^2].
       *}
      fGM2: extended;
      {*
       * Mass parameter for Earth-Moon barycenter GM in [AU^3/day^2].
       *}
      fGMB: extended;
      {*
       * Mass parameter for Mars GM in [AU^3/day^2].
       *}
      fGM4: extended;
      {*
       * Mass parameter for Jupiter GM in [AU^3/day^2].
       *}
      fGM5: extended;
      {*
       * Mass parameter for Saturn GM in [AU^3/day^2].
       *}
      fGM6: extended;
      {*
       * Mass parameter for Uranus GM in [AU^3/day^2].
       *}
      fGM7: extended;
      {*
       * Mass parameter for Neptune GM in [AU^3/day^2].
       *}
      fGM8: extended;
      {*
       * Mass parameter for Pluto GM in [AU^3/day^2].
       *}
      fGM9: extended;
      {*
       * Mass parameter for Sun GM in [AU^3/day^2].
       *}
      fGMS: extended;

      { GROUP 1050 from file header.xxx }
      {*
       * For each planet (and the Moon makes 10, and the Sun makes 11),Earth
       * nutations and Moon librations, each interval contains several complete
       * sets of coefficients, each covering a fraction of the interval duration.
       *}
      fNumberOfCoefSets: TDEIntArr;
      {*
       * Each planet (and the Moon makes 10, and the Sun makes 11), Earth
       * nutations and Moon librations, has a different number of Chebyshev
       * coefficients used to calculate each component of position and velocity.
       *}
      fNumberOfCoefs: TDEIntArr;
      {*
       * Initialize array for number of Chebyshev polynomials.
       *}
      fNumberOfPoly: TDEIntArr;
      {*
       * Define array for Chebyshev coefficients in one file. Real size will be
       * compute in subclasses.
       *}
      fEphemerisCoefficients: TDEExtArr;
      {*
       * Define array for start and final dates of file in Julian Days.
       *}
      fEphemerisDates: TDEExtArr;
      {*
       * Path to ASCII ephemerides files.
       *}
      fPathEph: String;
      procedure setPathEph(value: String);
    public
      property numbersPerInterval: Integer read fNumbersPerInterval;
      property startEpoch: Extended read fStartEpoch;
      property finalEpoch: Extended read fFinalEpoch;
      property intervalDuration: Integer read fIntervalDuration;
      property CLight: Extended read fCLight;
      property AU: Extended read fAU;
      property EMrat: Extended read fEMrat;
      property GM1: Extended read fGM1;
      property GM2: Extended read fGM2;
      property GMB: Extended read fGMB;
      property GM4: Extended read fGM4;
      property GM5: Extended read fGM5;
      property GM6: Extended read fGM6;
      property GM7: Extended read fGM7;
      property GM8: Extended read fGM8;
      property GM9: Extended read fGM9;
      property GMS: Extended read fGMS;
      property numberOfCoefSets: TDEIntArr read fNumberOfCoefSets;
      property numberOfCoefs: TDEIntArr read fNumberOfCoefs;
      property numberOfPoly: TDEIntArr read fNumberOfPoly;
      property ephemerisCoefficients: TDEExtArr read fEphemerisCoefficients;
      property ephemerisDates: TDEExtArr read fEphemerisDates;
      property pathEph: String read fPathEph write setPathEph;

      constructor Create;
      {*
       * @return string with planetary ephemeris version.
       *}
      function toString: String; override;
      {*
       * Prototype for the methods to read files with Chebyshev coefficients.
       * Method is override in the subclasses for the specific ephemeris.
       *
       * @param julTime
       *            Julian date for calculation.
       * @return true if the operations of reading is success.
       *}
      function readEphCoeff(julTime: Extended): Boolean; virtual; abstract;
  end;

implementation

constructor TDEheader.Create;
const
  _NumberOfCoefSets: array[0..13] of Integer = (0, 4, 2, 2, 1, 1, 1, 1, 1, 1, 8, 2, 4, 4);
  _NumberOfCoefs: array[0..13] of Integer = (0, 14, 10, 13, 11, 8, 7, 6, 6, 6, 13, 11, 10, 10);
  _NumberOfPoly: array[0..13] of Integer = (0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3);
begin
  // fNumbersPerInterval:=1016;
  fNumberOfCoefSets := _NumberOfCoefSets;
  fNumberOfCoefs := _NumberOfCoefs;
  fNumberOfPoly := _NumberOfPoly;
  SetLength(fEphemerisDates, 3);
  fPathEph := '';
end;

procedure TDEheader.setPathEph(value: String);
begin
  fPathEph := value;
end;

function TDEheader.toString: String;
begin
  Result := 'JPL Planetary Ephemeris DE'+IntToStr(fDEnomber)+'/LE'+IntToStr(fDEnomber);
end;

end.

