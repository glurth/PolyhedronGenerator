using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class SpinOnAxis : MonoBehaviour
{
    public Vector3 axis;
    public float degPerSec;

    void Update()
    {
        Quaternion change =Quaternion.AngleAxis(degPerSec * Time.deltaTime, axis);
        transform.localRotation *= change;
    }
}
