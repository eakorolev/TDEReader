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
      class function getPlanetPosvel(de: TDEheader; julTime: Double; i: Integer): TDEExtArr;
  end;

implementation

class function TDEephem.getPlanetPosvel(de: TDEheader; julTime: Double; i: Integer): TDEExtArr;
var
  interval, numbersToSkip, pointer, j, k, subInterval: Integer;
  intervalStartTime, subIntervalDuration, chebyshevTime: Extended;
  ephemerisRV, positionPoly, velocityPoly: TDEExtArr;
  coef: TDEExtMatr;
begin
  SetLength(ephemerisRV, 7);

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
      ephemerisRV[j] := 0;
      for k := 1 to de.numberOfCoefs[i] do
        ephemerisRV[j] := ephemerisRV[j] + coef[j][k] * positionPoly[k];

      {* Convert from km to A.U. *}
      if i < 12 then
        ephemerisRV[j] := ephemerisRV[j] / de.AU;
    end;

    {* Calculate the Chebyshev velocity polynomials. *}
    velocityPoly[1] := 0;
    velocityPoly[2] := 1;
    velocityPoly[3] := 4 * chebyshevTime;
    for j := 4 to de.numberOfCoefs[i] do
      velocityPoly[j] := 2 * chebyshevTime * velocityPoly[j - 1] + 2 * positionPoly[j - 1] - velocityPoly[j - 2];

    {* Calculate the velocity of the i'th planet. *}
    for j := de.numberOfPoly[i] + 1 to 2 * de.numberOfPoly[i] do begin
      ephemerisRV[j] := 0;
      for k := 1 to de.numberOfCoefs[i] do
        ephemerisRV[j] := ephemerisRV[j] + coef[j - de.numberOfPoly[i]][k] * velocityPoly[k];
      {*
       * The next line accounts for differentiation of the iterative
       * formula with respect to chebyshev time. Essentially, if dx/dt
       * = (dx/dct) times (dct/dt), the next line includes the factor
       * (dct/dt) so that the units are km/day.
       *}
      ephemerisRV[j] := ephemerisRV[j] * (2.0 * de.numberOfCoefSets[i] / de.intervalDuration);

      {* Convert from km to A.U. *}
      if i < 12 then
        ephemerisRV[j] := ephemerisRV[j] / de.AU;
    end;
  end else begin
    for j := 1 to 6 do
      ephemerisRV[j] := NaN;
  end;

  Result := ephemerisRV;
end;

end.

