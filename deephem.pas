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
     * @param DE
     *            class containing constants and the coefficients of Chebyshev
     *            polynomials for a specific JPL ephemeris.
     * @param JulTime
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
    class function GetPlanetPosvel(DE: TDEheader; JulTime: extended;
      i: integer): TDEExtArr;

    {*
     * Method gives the position and velocity of the point 'ntarg' with respect
     * to 'ncent'. This method is Java version of the procedure 'PLEPH', which
     * not repeats fortran source exactly, but return a similar result.
     *
     * @param DE
     *            Class containing the parameters of the JPL ephemeris and
     *            coefficients of the Chebyshev polynomials.
     * @param JulTime
     *            julian ephemeris date at which interpolation is wanted.
     * @param nTarg
     *            integer number of target point: 1 = MERCURY, 2 = VENUS, 3 =
     *            EARTH, 4 = MARS, 5 = JUPITER, 6 = SATURN , 7 = URANUS, 8 =
     *            NEPTUNE, 9 = PLUTO, 10 = MOON, 11 = SUN, 12 = SOLAR-SYSTEM
     *            BARYCENTER, 13 = EARTH-MOON BARYCENTER, 14 = NUTATIONS
     *            (LONGITUDE AND OBLIQ), 15 = LIBRATIONS, IF ON EPH FILE.
     * @param nCent
     *            integer number of center point, same as 'ntarg'.
     * @return array with posicions in [AU] and velocity in [AU/day]. In the
     *         case of nutations the first four number will be set to nutations
     *         and rate in [rad] and [rad/day] respectively. If 'ntarg=15' in
     *         [rad] and [rad/day].
     *}
    class function PlEph(DE: TDEheader; JulTime: extended;
      nTarg, nCent: integer): TDEExtArr;
  end;

implementation

class function TDEephem.GetPlanetPosvel(DE: TDEheader; JulTime: extended;
  i: integer): TDEExtArr;
var
  Interval, NumbersToSkip, Pointer, j, k, SubInterval: integer;
  IntervalStartTime, SubIntervalDuration, ChebyshevTime: extended;
  PositionPoly, VelocityPoly: TDEExtArr;
  Coef: TDEExtMatr;
begin
  SetLength(Result, 7);

  SetLength(PositionPoly, 20);
  SetLength(Coef, 4, 20);
  SetLength(VelocityPoly, 20);

  if (DE.ReadEphCoeff(JulTime)) and (i <= 13) and (i > 0) then
  begin
    Interval := floor((JulTime - DE.EphemerisDates[1]) / DE.IntervalDuration) + 1;
    IntervalStartTime := (Interval - 1) * DE.IntervalDuration + DE.EphemerisDates[1];
    if DE.NumberOfCoefSets[i] <> 0 then
      SubIntervalDuration := DE.IntervalDuration / DE.NumberOfCoefSets[i];
    SubInterval := floor((JulTime - IntervalStartTime) / SubIntervalDuration) + 1;
    NumbersToSkip := (Interval - 1) * DE.NumbersPerInterval;

    {*
     * Starting at the beginning of the coefficient array, skip the
     * first "numbers_to_skip" coefficients. This puts the pointer on
     * the first piece of data in the correct interval.
     *}
    Pointer := NumbersToSkip + 1;

    // Skip the coefficients for the first (i-1) planets
    for j := 1 to (i - 1) do
      Pointer := Pointer + DE.NumberOfPoly[j] * DE.NumberOfCoefSets[j] *
        DE.NumberOfCoefs[j];

    // Skip the next (subinterval - 1)*number_of_poly(i)*number_of_coefs(i) coefficients
    Pointer := Pointer + (SubInterval - 1) * DE.NumberOfPoly[i] * DE.NumberOfCoefs[i];

    for j := 1 to DE.NumberOfPoly[i] do
    begin
      for k := 1 to DE.NumberOfCoefs[i] do
      begin
        // Read the pointer'th coefficient as the array entry coef[j][k]
        Coef[j][k] := DE.EphemerisCoefficients[Pointer];
        Pointer := Pointer + 1;
      end;
    end;

    // Calculate the chebyshev time within the subinterval, between -1 and +1.
    ChebyshevTime := 2 * (JulTime - ((SubInterval - 1) * SubIntervalDuration +
      IntervalStartTime)) / SubIntervalDuration - 1;

    // Calculate the Chebyshev position polynomials.
    PositionPoly[1] := 1;
    PositionPoly[2] := ChebyshevTime;
    for j := 3 to DE.NumberOfCoefs[i] do
      PositionPoly[j] := 2 * ChebyshevTime * PositionPoly[j - 1] - PositionPoly[j - 2];

    // Calculate the position of the i'th planet at jultime.
    for j := 1 to DE.NumberOfPoly[i] do
    begin
      Result[j] := 0;
      for k := 1 to DE.NumberOfCoefs[i] do
        Result[j] := Result[j] + Coef[j][k] * PositionPoly[k];

      // Convert from km to A.U.
      if i < 12 then
        Result[j] := Result[j] / DE.AU;
    end;

    // Calculate the Chebyshev velocity polynomials.
    VelocityPoly[1] := 0;
    VelocityPoly[2] := 1;
    VelocityPoly[3] := 4 * ChebyshevTime;
    for j := 4 to DE.NumberOfCoefs[i] do
      VelocityPoly[j] := 2 * ChebyshevTime * VelocityPoly[j - 1] +
        2 * PositionPoly[j - 1] - VelocityPoly[j - 2];

    // Calculate the velocity of the i'th planet.
    for j := DE.NumberOfPoly[i] + 1 to 2 * DE.NumberOfPoly[i] do
    begin
      Result[j] := 0;
      for k := 1 to DE.NumberOfCoefs[i] do
        Result[j] := Result[j] + Coef[j - DE.NumberOfPoly[i]][k] * VelocityPoly[k];
      {*
       * The next line accounts for differentiation of the iterative
       * formula with respect to chebyshev time. Essentially, if dx/dt
       * = (dx/dct) times (dct/dt), the next line includes the factor
       * (dct/dt) so that the units are km/day.
       *}
      Result[j] := Result[j] * (2.0 * DE.NumberOfCoefSets[i] / de.IntervalDuration);

      // Convert from km to A.U.
      if i < 12 then
        Result[j] := Result[j] / DE.AU;
    end;
  end
  else
  begin
    for j := 1 to 6 do
      Result[j] := NaN;
  end;
end;

class function TDEephem.PlEph(DE: TDEheader; JulTime: extended;
  nTarg, nCent: integer): TDEExtArr;
var
  Targ, Cent, Earth, Moon: TDEExtArr;
  i: integer;
begin
  SetLength(Targ, 7);
  case nTarg of
    3:
    begin
      Earth := TDEephem.GetPlanetPosvel(DE, JulTime, 3);
      Moon := TDEephem.GetPlanetPosvel(DE, JulTime, 10);
      for i := 1 to 6 do
        Targ[i] := Earth[i] - Moon[i] / (1.0 + DE.EMrat);
    end;
    10:
    begin
      Earth := TDEephem.GetPlanetPosvel(DE, JulTime, 3);
      Moon := TDEephem.GetPlanetPosvel(DE, JulTime, 10);
      for i := 1 to 6 do
        Targ[i] := Earth[i] + Moon[i] - Moon[i] / (1.0 + DE.EMrat);
    end;
    12:
    begin
    end;
    13: Targ := TDEephem.GetPlanetPosvel(DE, JulTime, 3);
    14: Targ := TDEephem.GetPlanetPosvel(DE, JulTime, 12);
    15: Targ := TDEephem.GetPlanetPosvel(DE, JulTime, 13);
    else
      Targ := TDEephem.GetPlanetPosvel(DE, JulTime, nTarg);
  end;

  SetLength(Cent, 7);
  case nCent of
    0, 12, 14, 15:
    begin
    end;
    3:
    begin
      Earth := TDEephem.GetPlanetPosvel(DE, JulTime, 3);
      Moon := TDEephem.GetPlanetPosvel(DE, JulTime, 10);
      for i := 1 to 6 do
        Cent[i] := Earth[i] - Moon[i] / (1.0 + DE.EMrat);
    end;
    10:
    begin
      Earth := TDEephem.GetPlanetPosvel(DE, JulTime, 3);
      Moon := TDEephem.GetPlanetPosvel(DE, JulTime, 10);
      for i := 1 to 6 do
        Cent[i] := Earth[i] + Moon[i] - Moon[i] / (1.0 + DE.EMrat);
    end;
    13: Cent := TDEephem.GetPlanetPosvel(DE, JulTime, 3);
    else
      Cent := TDEephem.GetPlanetPosvel(DE, JulTime, nCent);
  end;

  SetLength(Result, 7);
  for i := 1 to 6 do
    Result[i] := Targ[i] - Cent[i];
end;

end.
