# nekdata
Data export and import to and from Nek5000/NekRS

Not sure where this project would go, but here are the TODO list:

1. Extra Nek5000/NekRS data into a box data structure.
   - Running Nek on a box mesh, interpolate the solution onto an unifom but non-overlappinh grid.
   - Read the outcome in python / matlab and fill the data into a `nx`-by-`ny` box
   - With the flexibility of incluing or excluding boundary points.
2. collection of nek reader / writer
