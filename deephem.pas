{* TDEReader: Peter Hristozov's JDEread (http://sourceforge.net/projects/jderead/) ported to Pascal.
 * 
 * Ported by Evgeny Korolev E-mail: mail@eakorolev.ru
 *}
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
unit DEephem;

interface

uses
  Classes, SysUtils, Math, DEheader;

type
  {*
   * Class contains a method for calculating the position and velocity of the
   * eight major planets and Pluto, Moon and Sun, also earth nutations and lunar
   * librations for each JPL ephemerides.
   *
   * @author Peter Hristozov E-mail: peterhri@hotmail.com
   * @version 1.3 2011 Jun 15.
   *}
  TDEephem = class
    public
      {*
       * Method calculate the position and velocity of body 'i' (nutations and
       * librations also) in epoch 'jultime' for a specific JPL ephemeris.<br>
       * The general idea is as follows: First read the coefficients corresponding
       * to 'jultime', and then calculate the positions and velocities of the
       * planet.
       *
       * @param de
       *            class containing constants and the coefficients of Chebyshev
       *            polynomials for a specific JPL ephemeris.
       * @param jultime
       *            Julian ephemeris epoch at which interpolation is wonted.
       * @param i
       *            designate of the astronomical bodies: 1 Mercury, 2 Venus, 3
       *            Earth-Moon barycenter, 4 Mars, 5 Jupiter, 6 Saturn, 7 Uranus,
       *            8 Neptune, 9 Pluto, 10 Geocentric Moon, 11 Sun, 12 nutations
       *            in longitude and obliquity (if on file), 13 lunar librations
       *            (if on file).
       * @return array with solar system barycentric positions and velocity for
       *         planets and Sun or geocentric position and velocity for Moon, for
       *         all in units [AU] and [AU/day]. The order of components is: X, Y,
       *         Z, dX, dY, dZ. If i=12 first four components return Earth
       *         nutation in longitude and obliquity and their rates in units
       *         [rad] and [rad/day]. If i=13 array contains lunar librations
       *         angles phi, thita and psi and their rates in units [rad] and
       *         [rad/day].
       *}
      class function GetPlanetPosvel(de: TDEheader; julTime: Extended; i: Integer): TDEExtArr;

      {*
       * Method gives the position and velocity of the point 'ntarg' with respect
       * to 'ncent'. This method is Java version of the procedure 'PLEPH', which
       * not repeats fortran source exactly, but return a similar result.
       *
       * @param de
       *            Class containing the parameters of the JPL ephemeris and
       *            coefficients of the Chebyshev polynomials.
       * @param et
       *            julian ephemeris date at which interpolation is wanted.
       * @param ntarg
       *            integer number of target point: 1 = MERCURY, 2 = VENUS, 3 =
       *            EARTH, 4 = MARS, 5 = JUPITER, 6 = SATURN , 7 = URANUS, 8 =
       *            NEPTUNE, 9 = PLUTO, 10 = MOON, 11 = SUN, 12 = SOLAR-SYSTEM
       *            BARYCENTER, 13 = EARTH-MOON BARYCENTER, 14 = NUTATIONS
       *            (LONGITUDE AND OBLIQ), 15 = LIBRATIONS, IF ON EPH FILE.
       * @param ncent
       *            integer number of center point, same as 'ntarg'.
       * @return array with posicions in [AU] and velocity in [AU/day]. In the
       *         case of nutations the first four number will be set to nutations
       *         and rate in [rad] and [rad/day] respectively. If 'ntarg=15' in
       *         [rad] and [rad/day].
       *}
      class function PlEph(de: TDEheader; et: Extended; nTarg, nCent: Integer): TDEExtArr;
  end;

implementation

class function TDEephem.GetPlanetPosvel(DE: TDEheader; JulTime: Extended; i: Integer): TDEExtArr;
var
  Interval, NumbersToSkip, Pointer, j, k, SubInterval: Integer;
  IntervalStartTime, SubIntervalDuration, ChebyshevTime: Extended;
  PositionPoly, VelocityPoly: TDEExtArr;
  Coef: TDEExtMatr;
begin
  SetLength(Result, 7);

  SetLength(positionPoly, 20);
  SetLength(coef, 4, 20);
  SetLength(velocityPoly, 20);

  if (de.readEphCoeff(julTime)) and (i <= 13) and (i > 0) then begin
    interval := floor((julTime - de.ephemerisDates[1]) / de.intervalDuration) + 1;
    intervalStartTime := (interval - 1) * de.intervalDuration + de.ephemerisDates[1];
    if de.numberOfCoefSets[i] <> 0 then
      subIntervalDuration := de.intervalDuration / de.numberOfCoefSets[i];
    subInterval := floor((julTime - intervalStartTime) / subIntervalDuration) + 1;
    numbersToSkip := (interval - 1) * de.numbersPerInterval;

    {*
     * Starting at the beginning of the coefficient array, skip the
     * first "numbers_to_skip" coefficients. This puts the pointer on
     * the first piece of data in the correct interval.
     *}
    pointer := numbersToSkip + 1;

    {* Skip the coefficients for the first (i-1) planets *}
    for j:=1 to (i - 1) do
      pointer := pointer + de.numberOfPoly[j] * de.numberOfCoefSets[j] * de.numberOfCoefs[j];

    {*
     * Skip the next (subinterval - 1)*number_of_poly(i)*number_of_coefs(i) coefficients
     *}
    pointer := pointer + (subInterval - 1) * de.numberOfPoly[i] * de.numberOfCoefs[i];

    for j := 1 to de.numberOfPoly[i] do begin
      for k := 1 to de.numberofcoefs[i] do begin
        {*
         * Read the pointer'th coefficient as the array entry
         * coef[j][k]
         *}
        coef[j][k] := de.ephemerisCoefficients[pointer];
        pointer := pointer + 1;
      end;
    end;

    {*
     * Calculate the chebyshev time within the subinterval, between -1
     * and +1.
     *}
    chebyshevTime := 2 * (julTime - ((subInterval - 1) * subIntervalDuration + intervalStartTime)) / subIntervalDuration - 1;

    {* Calculate the Chebyshev position polynomials. *}
    positionPoly[1] := 1;
    positionPoly[2] := chebyshevTime;
    for j := 3 to de.numberofcoefs[i] do
      positionPoly[j] := 2 * chebyshevTime * positionPoly[j - 1] - positionPoly[j - 2];

    {* Calculate the position of the i'th planet at jultime. *}
    for j := 1 to de.numberOfPoly[i] do begin
      Result[j] := 0;
      for k := 1 to de.numberOfCoefs[i] do
        Result[j] := Result[j] + coef[j][k] * positionPoly[k];

      {* Convert from km to A.U. *}
      if i < 12 then
        Result[j] := Result[j] / de.AU;
    end;

    {* Calculate the Chebyshev velocity polynomials. *}
    velocityPoly[1] := 0;
    velocityPoly[2] := 1;
    velocityPoly[3] := 4 * chebyshevTime;
    for j := 4 to de.numberOfCoefs[i] do
      velocityPoly[j] := 2 * chebyshevTime * velocityPoly[j - 1] + 2 * positionPoly[j - 1] - velocityPoly[j - 2];

    {* Calculate the velocity of the i'th planet. *}
    for j := de.numberOfPoly[i] + 1 to 2 * de.numberOfPoly[i] do begin
      Result[j] := 0;
      for k := 1 to de.numberOfCoefs[i] do
        Result[j] := Result[j] + coef[j - de.numberOfPoly[i]][k] * velocityPoly[k];
      {*
       * The next line accounts for differentiation of the iterative
       * formula with respect to chebyshev time. Essentially, if dx/dt
       * = (dx/dct) times (dct/dt), the next line includes the factor
       * (dct/dt) so that the units are km/day.
       *}
      Result[j] := Result[j] * (2.0 * de.numberOfCoefSets[i] / de.intervalDuration);

      {* Convert from km to A.U. *}
      if i < 12 then
        Result[j] := Result[j] / de.AU;
    end;
  end else begin
    for j := 1 to 6 do
      Result[j] := NaN;
  end;
end;

class function TDEephem.PlEph(de: TDEheader; et: Extended; nTarg, nCent: Integer): TDEExtArr;
var
  targ, cent, earth, moon: TDEExtArr;
  i: Integer;
begin
  SetLength(targ, 7);
  case nTarg of
    3: begin
      earth := TDEephem.GetPlanetPosvel(de, et, 3);
      moon := TDEephem.GetPlanetPosvel(de, et, 10);
      for i:= 1 to 6 do targ[i] := earth[i] - moon[i] / (1.0 + de.EMrat);
    end;
    10: begin
      earth := TDEephem.GetPlanetPosvel(de, et, 3);
      moon := TDEephem.GetPlanetPosvel(de, et, 10);
      for i:= 1 to 6 do targ[i] := earth[i] + moon[i] - moon[i] / (1.0 + de.EMrat);
    end;
    12: begin end;
    13: targ := TDEephem.GetPlanetPosvel(de, et, 3);
    14: targ := TDEephem.GetPlanetPosvel(de, et, 12);
    15: targ := TDEephem.GetPlanetPosvel(de, et, 13);
    else targ := TDEephem.GetPlanetPosvel(de, et, nTarg);
  end;

  SetLength(cent, 7);
  case nCent of
    0: begin end;
    3: begin
      earth := TDEephem.GetPlanetPosvel(de, et, 3);
      moon := TDEephem.GetPlanetPosvel(de, et, 10);
      for i:= 1 to 6 do cent[i] := earth[i] - moon[i] / (1.0 + de.EMrat);
    end;
    10: begin
      earth := TDEephem.GetPlanetPosvel(de, et, 3);
      moon := TDEephem.GetPlanetPosvel(de, et, 10);
      for i:= 1 to 6 do cent[i] := earth[i] + moon[i] - moon[i] / (1.0 + de.EMrat);
    end;
    12: begin end;
    13: cent := TDEephem.GetPlanetPosvel(de, et, 3);
    14: begin end;
    15: begin end;
    else cent := TDEephem.GetPlanetPosvel(de, et, nCent);
  end;

  SetLength(Result, 7);
  for i:= 1 to 6 do Result[i] := targ[i] - cent[i];
end;

end.

